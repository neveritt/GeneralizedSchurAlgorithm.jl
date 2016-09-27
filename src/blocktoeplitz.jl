import ToeplitzMatrices: AbstractToeplitz, Toeplitz
import Base: size, getindex, convert
"""
A block Toeplitz matrix is constructed by the first block column and
first block row.

# Examples
```julia
julia>T0 = randn(2,3)
vc = randn(10,3)
vr = randn(2,8)
vc[1:2,1:3] = T0
vr[1:2,1:3] = T0
T = Toeplitz(vc,vr)
```
"""
#= Block Toeplitz matrix
 A Block Toeplitz matrix T ∈ R^{mk,nl} is a matrix that has blocks T_i of size
 k by l collected in a Toeplitz structur:
  [T_0       T_1   …   T\_{n-1}
   T_{-1}   ⋱         ⋮
   T_{-2}   ⋱   ⋱    ⋮
   T_{-m+1} …    …   T_{-m+n}]=#

# General BlockToeplitz matrix
immutable BlockToeplitz{T<:AbstractFloat} <: AbstractToeplitz{T}
    vc::Matrix{T}
    vr::Matrix{T}

    @compat function (::Type{BlockToeplitz}){T<:AbstractFloat}(vc::Matrix{T}, vr::Matrix{T})
      k = size(vr,1)
      l = size(vc,2)
      if !isapprox(vc[1:k,1:l], vr[1:k,1:l])
        warn("First block element must be the same")
        throw(DomainError())
      end
      new{T}(vc,vr)
    end
end

# constructor
Toeplitz(col::AbstractMatrix, row::AbstractMatrix) = BlockToeplitz(col, row)

# Size of a general Toeplitz matrix
function size(A::BlockToeplitz, dim::Int)
  if dim == 1
    return size(A.vc,1)
  elseif dim == 2
    return size(A.vr,2)
  elseif dim > 2
    return 1
  else
    warn("arraysize: dimension out of range")
    throw(DomainError())
  end
end

# Blocksize of a general Toeplitz matrix
function blocksize(A::BlockToeplitz, dim::Int)
  if dim == 1
    return convert(Int,div(size(A.vc,1),sizeofblock(A,1)))
  elseif dim == 2
    return convert(Int,div(size(A.vr,2),sizeofblock(A,2)))
  elseif dim > 2
    return 1
  else
    warn("arraysize: dimension out of range")
    throw(DomainError())
  end
end
blocksize(A::BlockToeplitz) = (blocksize(A,1), blocksize(A,2))

# Blocksize of a general Toeplitz matrix
function sizeofblock(A::BlockToeplitz, dim::Int)
  if dim == 1
    return size(A.vr,1)
  elseif dim == 2
    return size(A.vc,2)
  elseif dim > 2
    return 1
  else
    warn("arraysize: dimension out of range")
    throw(DomainError())
  end
end
sizeofblock(A::BlockToeplitz) = (sizeofblock(A,1), sizeofblock(A,2))

function getindex(A::BlockToeplitz, i::Integer, j::Integer)
  k,l = sizeofblock(A)
  m,n = blocksize(A)
  if i > size(A,1) || j > size(A,2)
    warn("BoundsError()")
    throw(DomainError())
  end
  blockidx = div(j-1,l) - div(i-1,k)
  if blockidx >= 0
    return A.vr[mod(i-1,k) + 1, blockidx*l + mod(j-1,l) + 1]
  else
    return A.vc[-blockidx*k + mod(i-1,k) + 1,  mod(j-1,l) + 1 ]
  end
end

function getblock(A::BlockToeplitz, i::Integer)
  m,n = blocksize(A)
  if i > n-1 || i < -m+1
    warn("BoundsError()")
    throw(DomainError())
  end
  k,l = sizeofblock(A)
  if i >= 0
    return A.vr[1:k, i*l+(1:l)]
  else
    return A.vc[-i*k+(1:k), 1:l]
  end
end

getcol(A::BlockToeplitz) = A.vc
getrow(A::BlockToeplitz) = A.vr

convert(::Type{Matrix}, A::BlockToeplitz) = full(A)

# Full version of a BlockToeplitz matrix
function full{T}(A::BlockToeplitz{T})
  m, n = size(A)
  Af = Array(T, m, n)
  for j = 1:n
    for i = 1:m
      Af[i,j] = A[i,j]
    end
  end
  return Af
end

# Application of a general Toeplitz matrix to a column vector
function A_mul_B!{T}(α::T, A::BlockToeplitz{T}, x::StridedVector{T}, β::T, y::StridedVector{T})
  m = size(A,1)
  n = size(A,2)
  if m != length(y)
    throw(DimensionMismatch(""))
  end
  if n != length(x)
    throw(DimensionMismatch(""))
  end
  y[:] *= β
  for j = 1:n
    tmp = α * x[j]
    for i = 1:m
      y[i] += tmp*A[i,j]
    end
  end
  return y
end

# Application of a general Toeplitz matrix to a general matrix
function A_mul_B!{T}(α::T, A::BlockToeplitz{T}, B::StridedMatrix{T}, β::T, C::StridedMatrix{T})
  l = size(B, 2)
  if size(C, 2) != l
      throw(DimensionMismatch("input and output matrices must have same number of columns"))
  end
  for j = 1:l
    A_mul_B!(α::T, A, view(B, :, j), β::T, view(C, :, j))
  end
  return C
end

(*){T}(A::BlockToeplitz{T}, B::StridedMatrix{T}) =
     A_mul_B!(one(T), A, B, zero(T), zeros(T, size(A,1), size(B,2)))

(*){T}(A::BlockToeplitz{T}, B::StridedVector{T}) =
      A_mul_B!(one(T), A, B, zero(T), zeros(T, size(A,1)))
