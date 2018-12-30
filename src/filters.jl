
# Univariate HP filter, Hodrick and Prescott, 1997, “Postwar U.S. Business Cycles: An Empirical Investigation.” Journal of Money, Credit, and Banking. Vol. 29, No. 1.

"""
    hp_filter(y::Array{T,1}, λ::Real)
    HP filter with parameter `λ`.

"""
function hpfilter(y::Array{Tp,1}, λ::Real) where Tp <: Real
    T = length(y)
    T < 3 ? error("Sample size should be at least 3") : nothing

    if T == 3
        A = I + λ*[1.0 -2.0  1.0
                  -2.0  4.0 -2.0
                   1.0 -2.0  1.0]
    else
        d0 = 1.0 .+ [λ; 5.0*λ; 6.0*λ*ones(T-4); 5.0*λ; λ]
        d1 = [-2.0*λ; -4.0*λ*ones(T-3); -2.0*λ]
        d2 = λ*ones(T-2)
         
        A  = spdiagm(-2 => d2, -1 => d1, 0 => d0, 1 => d1, 2 => d2)
    end

    trend    = A\y
    cyclical = y - trend
    
    return cyclical,trend
end

##################################################################################
# Optimal HP Filter Parameter (Dermoune et al. 2008), Theorem 1
function optimalλDermoune(y)
    # Constructing the matrix P
    T = length(y)
    # P = zeros(n-2,n)
    # a = [1.0, -2.0, 1.0] 
    # @inbounds for i in 1:n-2
    #             P[i,i:i+2] = a
    #           end

    # py      = P*y
    py    = diff(diff(y)) 
    num   = dot(py,py)
    denum = dot(py[2:T-2],py[1:T-3])
    # The consistent estimator for λ          
    λopt = (-1/4)/(3/2 + ((T-3)*num)/((T-2)*denum))^-1
  
end

##################################################################################
# Pedersen 2001

function optimalλPedersen(y,wH)
    #w1      = pi/10    # cut off frequency based on average duration of business cycle
    # Optimal λ minimizes the distortion
    obj = x->distortion(y,wH,x)
    result = optimize(obj, 1.0,10000.0) # seaching on [1.0,10_000.0], which should be ok.

    if Optim.converged(result)
        return result.minimizer[1]
    else
        error("Couldn't find the minimizer...")
    end

end

function distortion(y,wH,λ)        
    T      = length(y)

    # normalized psd
    xdft   = fft(y)
    xdft   = xdft[1:Int(floor(T/2+1))]
    F      = (1/(2*pi*T)) * abs2.(xdft)
    F[2:end-1] = 2*F[2:end-1]
    w      = range(0.0,stop=pi,length=Int(floor(T/2+1)))

    dw      = w[2] - w[1]
    HIdeal  = zeros(length(w))
    HIdeal[abs.(w) .>= wH] .= 1.0
    HHP     = @. abs2((4.0*λ*(1.0 - cos.(w))^2)/(4.0*λ*(1.0 - cos(w))^2 + 1.0))
    Norm    = sum(2.0*F*dw)
    v       = 2.0*F*dw/Norm

    Q       = dot(abs.(HIdeal - HHP),v)
end
####################################################################################################################################################################
####################################################################################################################################################################

# Baxter-King filter, Baxter and King (1999) 

# some recomended values:
#--------------------------------------------------------
# Frequency	    Lead/lag	Lower limit	    Upper limit
# Annual	        3	        2	            8
# Semi-annual	    6	        3	            16
# Quarterly	        12	        6	            32
# Monthly	        36	        8	            96
# Weekly	        156	        78	            416
#--------------------------------------------------------

function bkfilter(y::Array{Tp,1},pL::Real,pU::Real,L::Int) where Tp <: Real

    pL<2.0 ? error("Lower cut-off cannot be smaller than 2!") : nothing
    pL>=pU ? error("Lower cut-off cannot be bigger than upper cut-off!") : nothing
    
    T = length(y)

    wL = 2.0*pi/pL
    wU = 2.0*pi/pU
    b0 = (wL - wU)/pi
    b  = zeros(L)

    @inbounds for j in 1:L 
        b[j] = (sin(wL*j) - sin(wU*j))/(j*pi)
    end

    meanB = (b0 + 2.0*sum(b))/(2*L + 1)

    b0 -= meanB
    b .-= meanB
    B = [reverse(b); b0; b]
    
    # will be losing L observation from start and end
    cyclical =  Array{Float64,1}(undef,T)
    cyclical[1:L]     .= NaN 
    cyclical[T+1-L:T] .= NaN

    @inbounds for t in (L + 1):(T - L)
        cyclical[t]  = dot(B,y[t-L:t+L])   
    end

    return cyclical

end
