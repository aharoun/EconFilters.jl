using Revise
using EconFilters
using BenchmarkTools

# random data
x = collect(1:1000) + randn(1000)

c,t  = hp_filter(x,1600)

λOpt       =  optimalλDermoune(x)
cOpt,tOpt  = hp_filter(x,λOpt)

