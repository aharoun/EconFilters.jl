
# HP filter, Hodrick and Prescott, 1997, “Postwar U.S. Business Cycles: An Empirical Investigation.” Journal of Money, Credit, and Banking. Vol. 29, No. 1.
function hp_filter(data::Array{T,1}, λ::Real) where T <: Real
    #naber
    n = length(data)
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

    trend    = A\data
    cyclical = data - trend
    
    return cyclical,trend
end

##################################################################################
# Optimal HP Filter Parameter (Dermoune et al. 2008), Theorem 1
function optimalλDermoune(data)
    # Constructing the matrix P
    n = length(data)
    P = zeros(n-2,n)
    a = [1.0, -2.0, 1.0] 
    @inbounds for i in 1:n-2
                P[i,i:i+2] = a
              end

    # The consistent estimator for λ
    py      = P*data
    num     = dot(py,py)
    denum   = dot(py[2:n-2],py[1:n-3])

    λopt = (-1/4)*(3/2 + ((n-3)*num)/((n-2)*denum)^-1)
  
end
##################################################################################

