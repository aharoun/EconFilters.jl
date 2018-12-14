using Revise
using EconFilters
using BenchmarkTools
using Plots
using CSV

# load US real GDP (quarterly)
data = CSV.read("example/GDPC1.csv")
y = log.(data[:GDPC1])

# HP with λ = 1600   
t,c  = hpfilter(y,1600)

# HP with optimal λ based on Dermoune et al (2008)
λOpt       =  optimalλDermoune(y)
tOpt,cOpt  = hpfilter(y,λOpt)

# HP with optimal λ based on Pedersen (2001)
λOpt2       =  optimalλPedersen(y,pi/10)  
tOpt2,cOpt2 = hpfilter(y,λOpt2)

plot(data[:date],[c cOpt cOpt2],label=["1600","$λOpt","$λOpt2"],legend=:bottomright,xtickfont = font(6))
plot(data[:date],[t tOpt tOpt2 y],label=["1600","$λOpt","$λOpt2","data"],legend=:bottomright,xtickfont = font(6))


# todo; one sided filter
# extended filter
# BK filter