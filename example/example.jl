using Revise
using EconFilters
using Plots
using CSV

# load US real GDP (quarterly)
data = CSV.read("example/GDPC1.csv")
y = log.(data[:GDPC1])

# HP with λ = 1600   
c,t  = hpfilter(y,1600)

# HP with optimal λ based on Pedersen (2001)
λOpt      =  optimalλPedersen(y,pi/10)  
cOpt,tOpt = hpfilter(y,λOpt)

# BK Filter 
# bkfilter(data,pL,pU,L) where pL and pU are lower and upper periods for bandpass filter, L is the order of the symmetric moving average
cBK =  bkfilter(y,6,32,12)

plot(data[:date],[c cOpt cBK],label=["Hp1600","Hp$λOpt","BK"],legend=:bottomright,xtickfont = font(5))

