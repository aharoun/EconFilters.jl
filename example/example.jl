using Revise
using EconFilters
using Plots
using CSV

# load US real GDP (quarterly)
data = CSV.read("GDPC1.csv")
y = log.(data[:GDPC1])

# HP with λ = 1600   
c,t  = hpfilter(y,1600)

# HP with optimal λ based on Pedersen (2001)
λOpt      =  optimalλPedersen(y,pi/10)  
cOpt,tOpt = hpfilter(y,λOpt)

# BK Filter 
# bkfilter(data,pL,pU,L) where pL and pU are lower and upper periods for bandpass filter, L is the order of the symmetric moving average
cBK =  bkfilter(y,6,32,12)


# To reduce end-of-sample bias of hp filter: Estimate an arima model, extend y series with forecasted values.
# For this, use SSM package (unregistered).
# First add the package
# ```add https://github.com/aharoun/SSM.jl```
#
using SSM
model = arima(2,1,2)
modelEst, estParams , _ = estimate(model, y);

fH = 10  # forecast horizon (period)
yF, _ = forecast(modelEst, y, fH)

cExt, tExt = hpfilter([y;yF],1600)

plot(data[:date],[c cExt[1:end-fH] cOpt cBK],label=["Hp1600","Hp1600 Extended","Hp$λOpt","BK"],legend=:bottomright,xtickfont = font(5))

