module EconFilters
using LinearAlgebra
using SparseArrays
using FFTW
using Optim



include("hpfilter.jl")
#include("baxterking_filter.jl")

export hpfilter,optimalλDermoune,distortion,optimalλPedersen



end