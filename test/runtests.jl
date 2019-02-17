using EconFilters, Test

y = collect(1:100) + randn(100)
# HP with λ = 1600
c,t  = hpfilter(y,1600)
@test any(.!isnan.(c))

# Optimal HP Filter Parameter (Dermoune et al. 2008)
λOpt1 =  optimalλDermoune(y)
# HP with optimal λ based on Pedersen (2001)
λOpt2 =  optimalλPedersen(y,pi/10)

# BK Filter
# bkfilter(data,pL,pU,L) where pL and pU are lower and upper periods for bandpass filter, L is the order of the symmetric moving average
cBK =  bkfilter(y,6,32,12)

