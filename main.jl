using Revise
using EconFilters
using BenchmarkTools
using Plots

# random data
n = 50
y = collect(range(1,stop=10.0,length = n)) + randn(n)

# 
c,t  = hp_filter(y,1600)

λOpt       =  optimalλDermoune(y)
cOpt,tOpt  = hp_filter(y,λOpt)

λOpt2       =  optimalλPedersen(y,pi/10)
cOpt2,tOpt2  = hp_filter(y,λOpt2)

plot(c)
plot!(cOpt)
plot!(cOpt2)

