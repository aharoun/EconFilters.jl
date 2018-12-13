module EconFilters
using LinearAlgebra
using SparseArrays



include("hp_filter.jl")
#include("baxterking_filter.jl")

export hp_filter,optimalλDermoune



end