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
immutable BlockToeplitz{T<:Number,M1<:AbstractMatrix,M2<:AbstractMatrix} <: AbstractToeplitz{T}
    vc::M1
    vr::M2
    m::Int
    n::Int

    @compat function (::Type{BlockToeplitz}){M1<:AbstractMatrix,M2<:AbstractMatrix}(vc::M1, vr::M2)
      vc,vr = promote(vc,vr)
      k = size(vr,1)
      l = size(vc,2)
      mk = size(vc,1)
      nl = size(vr,2)
      m = convert(Int, round(mk/k))
      n = convert(Int, round(nl/l))
      if !isapprox(vc[1:k,1:l], vr[1:k,1:l])
        warn("First block element must be the same")
        throw(DomainError())
      end
      new{eltype(vc),typeof(vc),typeof(vr)}(vc,vr,m,n)
    end
end

# constructor
Toeplitz(col::AbstractMatrix, row::AbstractMatrix) = BlockToeplitz(col, row)

# Sizes a general Toeplitz matrix
function size(A::BlockToeplitz, dim::Int)
  @assert dim > 0 "size: dim must be positive"
  dim == 1 ? size(A.vc,1) : dim == 2 ? size(A.vr,2) : 1
end
size(A::BlockToeplitz)               = (size(A.vc,1), size(A.vr,2))
size(A::BlockToeplitz, dims::Int...) = map(x-> size(A, x), dims)

zero{T}(A::BlockToeplitz{T})         = BlockToeplitz(spzeros(T,size(A.vc)...),
  spzeros(T,size(A.vr)...))


eltype{T}(A::BlockToeplitz{T}) = T
promote_type{T,S}(A::BlockToeplitz{T}, c::S) = BlockToeplitz{promote_type(T,S)}

convert{T,S<:Number}(::Type{BlockToeplitz{T}}, c::S) = BlockToeplitz(fill(convert(T,c),1,1),fill(convert(T,c),1,1))

function blocksize(A::BlockToeplitz, dim::Int)
  @assert dim > 0 "blocksize: dim must be positive"
  dim == 1 ? A.m : dim == 2 ? A.n : 1
end
blocksize(A::BlockToeplitz)                 = (A.m, A.n)
blocksize(A::BlockToeplitz, dims::Int...)   = map(x-> blocksize(A, x), dims)

function sizeofblock(A::BlockToeplitz, dim::Int)
  @assert dim > 0 "sizeofblock: dim must be positive"
  dim == 1 ? size(A.vr,1) : dim == 2 ? size(A.vc,2) : 1
end
sizeofblock(A::BlockToeplitz)               = (size(A.vr,1), size(A.vc,2))
sizeofblock(A::BlockToeplitz, dims::Int...) = map(x-> sizeofblock(A, x), dims)

@compat Base.IndexStyle(::Type{<:BlockToeplitz}) = IndexCartesian()
getindex(A::BlockToeplitz, i::Int) = A[rem(i-1,size(A,1))+1, div(i-1,size(A,1))+1]

function checkbounds(A::BlockToeplitz, i::Integer, j::Integer)
  (i ≤ size(A,1) && j ≤ size(A,2)) || throw(BoundsError())
end

function getindex(A::BlockToeplitz, i::Int, j::Int)
  checkbounds(A, i, j)
  k,l = sizeofblock(A)
  blockidx = div(j-1,l) - div(i-1,k)
  if blockidx >= 0
    return A.vr[mod(i-1,k) + 1, blockidx*l + mod(j-1,l) + 1]
  end
  return A.vc[-blockidx*k + mod(i-1,k) + 1,  mod(j-1,l) + 1 ]
end

function checkblockbounds(A::BlockToeplitz, i::Integer)
  m,n = blocksize(A)
  (i > n-1 || i < -m+1) && throw(BoundsError())
end

function getblock(A::BlockToeplitz, i::Int)
  checkblockbounds(A, i)
  k,l = sizeofblock(A)
  if i >= 0
    @inbounds return A.vr[1:k, i*l+(1:l)]
  end
  @inbounds return A.vc[-i*k+(1:k), 1:l]
end

getcol(A::BlockToeplitz) = A.vc
getrow(A::BlockToeplitz) = A.vr

convert(::Type{Matrix}, A::BlockToeplitz) = full(A)

transpose(A::BlockToeplitz)  = BlockToeplitz(A.vr.', A.vc.')
ctranspose(A::BlockToeplitz) = BlockToeplitz(A.vr', A.vc')

# Full version of a BlockToeplitz matrix
function full{T}(A::BlockToeplitz{T})
  m, n = size(A)
  Af = Matrix{T}(m, n)
  @simd for i = 1:m
    @simd for j = 1:n
      @inbounds Af[i,j] = A[i,j]
    end
  end
  return Af
end

# Application of a general Toeplitz matrix to a column vector
function A_mul_B!{T}(α::T, A::BlockToeplitz{T}, x::StridedVector{T}, β::T, y::StridedVector{T})
  m, n = size(A)
  m == length(y) || throw(DimensionMismatch(""))
  n == length(x) || throw(DimensionMismatch(""))
  y[:] *= β
  @simd for j = 1:n
    @inbounds tmp = α * x[j]
    @simd for i = 1:m
      @inbounds y[i] += tmp*A[i,j]
    end
  end
  return y
end

# Application of a general Toeplitz matrix to a general matrix
function A_mul_B!{T}(α::T, A::BlockToeplitz{T}, B::StridedMatrix{T}, β::T, C::StridedMatrix{T})
  size(C, 2) == size(B, 2) || throw(DimensionMismatch("input and output matrices must have same number of columns"))
  @simd for j = 1:size(B, 2)
    @inbounds A_mul_B!(α::T, A, view(B, :, j), β::T, view(C, :, j))
  end
  return C
end

function A_mul_B_block!{T}(α::T, A::BlockToeplitz{T}, B::StridedMatrix{T}, β::T, C::StridedMatrix{T})
  size(C, 2) == size(B, 2) || throw(DimensionMismatch("input and output matrices must have same number of columns"))
  size(A, 2) == size(B, 1) || throw(DimensionMismatch("input and output matrices must have same number of columns"))
  k,l = sizeofblock(A)
  m,n = blocksize(A)
  C[:] .*= β
  @simd for i in 0:n-1
    Tᵢ = α.*getblock(A,i)
    @simd for col_idx in 1+i:min(n, m+i)
      row_idx   = col_idx-i
      C_row_idx = (row_idx-1)*k+(1:k)
      B_row_idx = (col_idx-1)*l+(1:l)
      @inbounds C[C_row_idx,:] += Tᵢ*view(B, B_row_idx, :)
    end
  end
  @simd for i in 1:m-1
    Tᵢ = α.*getblock(A,-i)
    @simd for col_idx in 1:min(m-i,n)
      row_idx   = col_idx+i
      C_row_idx = (row_idx-1)*k+(1:k)
      B_row_idx = (col_idx-1)*l+(1:l)
      @inbounds C[C_row_idx,:] += Tᵢ*view(B, B_row_idx, :)
    end
  end
  return C
end

+{T<:Number}(A1::BlockToeplitz{T}, A2::BlockToeplitz{T}) = BlockToeplitz(A1.vc+A2.vc, A1.vr+A2.vr)
-{T<:Number}(A1::BlockToeplitz{T}, A2::BlockToeplitz{T}) = BlockToeplitz(A1.vc-A2.vc, A1.vr-A2.vr)

(*){T}(A::BlockToeplitz{T}, B::StridedMatrix{T}) =
    A_mul_B_block!(one(T), A, B, zero(T), zeros(T, size(A,1), size(B,2)))

(*){T}(A::BlockToeplitz{T}, B::StridedVector{T}) =
    A_mul_B!(one(T), A, B, zero(T), zeros(T, size(A,1)))

# Application of a general Toeplitz matrix to a column vector
function At_mul_B!{T}(α::T, A::BlockToeplitz{T}, x::StridedVector{T}, β::T, y::StridedVector{T})
  n, m = size(A)
  if m != length(y)
    throw(DimensionMismatch(""))
  end
  if n != length(x)
    throw(DimensionMismatch(""))
  end
  y[:] *= β
  @simd for j = 1:n
    @inbounds tmp = α * x[j]
    @simd for i = 1:m
      @inbounds y[i] += tmp*A[j,i]
    end
  end
  return y
end

function At_mul_B_block!{T}(α::T, A::BlockToeplitz{T}, B::StridedMatrix{T}, β::T, C::StridedMatrix{T})
  size(C, 2) == size(B, 2) || throw(DimensionMismatch("input and output matrices must have same number of columns"))
  size(A, 1) == size(B, 1) || throw(DimensionMismatch(""))
  k,l = sizeofblock(A)
  m,n = blocksize(A)
  C[:] .*= β
  @simd for i in 0:m-1
    Tᵢ = α.*getblock(A,-i)
    @simd for col_idx in 1+i:min(m, n+i)
      row_idx   = col_idx-i
      C_row_idx = (row_idx-1)*l+(1:l)
      B_row_idx = (col_idx-1)*k+(1:k)
      @inbounds C[C_row_idx,:] += Tᵢ.'*view(B, B_row_idx, :)
    end
  end
  @simd for i in 1:n-1
    Tᵢ = α.*getblock(A,i)
    @simd for col_idx in 1:min(n-i,m)
      row_idx   = col_idx+i
      C_row_idx = (row_idx-1)*l+(1:l)
      B_row_idx = (col_idx-1)*k+(1:k)
      @inbounds C[C_row_idx,:] += Tᵢ.'*view(B, B_row_idx, :)
    end
  end
  return C
end

# Application of a general Toeplitz matrix to a general matrix
function At_mul_B!{T}(α::T, A::BlockToeplitz{T}, B::StridedMatrix{T}, β::T, C::StridedMatrix{T})
 l = size(B, 2)
 if size(C, 2) != l
     throw(DimensionMismatch("input and output matrices must have same number of columns"))
 end
  for j = 1:l
    @inbounds At_mul_B!(α::T, A, view(B, :, j), β::T, view(C, :, j))
  end
  return C
end

At_mul_B{T}(A::BlockToeplitz{T}, B::StridedMatrix{T}) =
           At_mul_B_block!(one(T), A, B, zero(T), zeros(T, size(A,2), size(B,2)))

At_mul_B{T}(A::BlockToeplitz{T}, B::StridedVector{T}) =
           At_mul_B!(one(T), A, B, zero(T), zeros(T, size(A,2)))
