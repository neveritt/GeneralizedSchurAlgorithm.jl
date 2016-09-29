module GeneralizedSchurAlgorithm

# Import functions for overloading
import Base: A_mul_B!, A_mul_Bc!
import Base: size, getindex, convert, ctranspose, (*)

import Base.LinAlg: BlasFloat, Char, BlasInt, LAPACKException, axpy!, BLAS.scal!

import ToeplitzMatrices: AbstractToeplitz, Toeplitz, full

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
  # h_rotation
  H_procedure,
  h_Algorithm,
  # oe_rotation
  OE_procedure,
  oe_algorithm,
  # generalized_schur
  schuralgorithm,
  # applications
  qrtoeplitz,
  lstoeplitz

# Polynomials package is needed
using Compat

# include files
include("blocktoeplitz.jl")
include("generalized_schur.jl")
include("h_rotation.jl")
include("oe_rotation.jl")
include("householder.jl")
include("lstoeplitz.jl")
include("qrtoeplitz.jl")

end # module
