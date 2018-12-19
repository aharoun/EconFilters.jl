
# Univariate HP filter, Hodrick and Prescott, 1997, “Postwar U.S. Business Cycles: An Empirical Investigation.” Journal of Money, Credit, and Banking. Vol. 29, No. 1.

"""
    hp_filter(y::Array{T,1}, λ::Real)
    HP filter with parameter `λ`.

"""
function hpfilter(y::Array{T,1}, λ::Real) where T <: Any
    n = length(y)
    n < 3 ? error("Sample size should be at least 3") : nothing

    if n == 3
        A = I + λ*[1.0 -2.0  1.0
                  -2.0  4.0 -2.0
                   1.0 -2.0  1.0]
    else
        d0 = 1.0 .+ [λ; 5.0*λ; 6.0*λ*ones(n-4); 5.0*λ; λ]
        d1 = [-2.0*λ; -4.0*λ*ones(n-3); -2.0*λ]
        d2 = λ*ones(n-2)
         
        A  = spdiagm(-2 => d2, -1 => d1, 0 => d0, 1 => d1, 2 => d2)
    end

    trend    = A\y
    cyclical = y - trend
    
    return trend,cyclical
end

##################################################################################
# Optimal HP Filter Parameter (Dermoune et al. 2008), Theorem 1
function optimalλDermoune(y)
    # Constructing the matrix P
    n = length(y)
    # P = zeros(n-2,n)
    # a = [1.0, -2.0, 1.0] 
    # @inbounds for i in 1:n-2
    #             P[i,i:i+2] = a
    #           end

    # py      = P*y
    py    = diff(diff(y)) 
    num   = dot(py,py)
    denum = dot(py[2:n-2],py[1:n-3])
    # The consistent estimator for λ          
    λopt = (-1/4)*(3/2 + ((n-3)*num)/((n-2)*denum))^-1
  
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
    n      = length(y)

    # normalized psd
    xdft   = fft(y)
    xdft   = xdft[1:Int(floor(n/2+1))]
    F      = (1/(2*pi*n)) * abs2.(xdft)
    F[2:end-1] = 2*F[2:end-1]
    w      = range(0.0,stop=pi,length=Int(floor(n/2+1)))

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