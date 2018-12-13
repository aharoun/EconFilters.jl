module EconFilters
using LinearAlgebra
using SparseArrays
using FFTW
using Optim



include("hp_filter.jl")
#include("baxterking_filter.jl")

export hp_filter,optimalλDermoune,distortion,optimalλPedersen



end