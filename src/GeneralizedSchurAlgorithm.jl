module GeneralizedSchurAlgorithm

# Import functions for overloading

import LinearAlgebra: mul!, rmul!, lmul!, norm, qr, triu, tril
import Base: size, eltype, getindex, convert, transpose, adjoint, promote_rule
import Compat.view

# this stuff still needs some more attention
import Base: promote_type, checkbounds
import LinearAlgebra:  Char, axpy!, BLAS.scal! #maybe not needed anymore
import LinearAlgebra: BLAS.@blasfunc, BlasFloat, BlasInt, LAPACKException, #maybe not everything here is needed
DimensionMismatch, SingularException, PosDefException, chkstride1, checksquare,
LAPACK.larf!, LAPACK.larfg!

import Base: one, zero, +, -, /, *
import ToeplitzMatrices: AbstractToeplitz, Toeplitz, Matrix

# Export only the useful functions
export
  # BlockToeplitz
  BlockToeplitz,
  Toeplitz,
  blocksize,
  sizeofblock,
  getblock,
  getcol,
  getrow,
  mul!,
  # h_rotation
  H_procedure,
  h_procedure,
  h_Algorithm,
  # oe_rotation
  OE_procedure,
  oe_procedure,
  oe_algorithm,
  # generalized_schur
  schuralgorithm,
  # applications
  qrtoeplitz,
  lstoeplitz

# Polynomials package is needed
using Compat

# include files
include("utilities.jl")
include("blocktoeplitz.jl")
include("generalized_schur.jl")
include("h_rotation.jl")
include("oe_rotation.jl")
include("householder.jl")
include("lstoeplitz.jl")
include("qrtoeplitz.jl")

end # module
