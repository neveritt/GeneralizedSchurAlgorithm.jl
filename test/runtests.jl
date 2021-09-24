using Test, LinearAlgebra, DSP, SparseArrays
import Compat.view
using GeneralizedSchurAlgorithm


# function to generate test-case
include("gen_toeplitz.jl")

# tests
include("blocktoeplitz.jl")
include("rotations.jl")
include("toeplitz.jl")
include("generalized_schur.jl")
