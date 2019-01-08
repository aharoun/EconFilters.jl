module EconFilters
using LinearAlgebra
using SparseArrays
using FFTW
using Optim



include("filters.jl")

export hpfilter,optimalλDermoune,optimalλPedersen,bkfilter



end
