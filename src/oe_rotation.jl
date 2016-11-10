

"""
    OE_procedure(i1,i2,α,β,ρ) -> H
A hyperbolic rotation operator implemented by the H-procedure.

The `OE_procedure` type supports left multiplication `H*A` and conjugated transpose
right multiplication `A*H'`. The type doesn't have a `size` and can therefore
be multiplied with matrices of arbitrary size as long as `i2<=size(A,2)` for
`H*A` or `i2<=size(A,1)` for `A*H`.
"""
immutable OE_procedure{T}
  i1::Int
  i2::Int
  c::T
  s::T
end

convert{T}(::Type{OE_procedure{T}}, H::OE_procedure{T}) = H
convert{T}(::Type{OE_procedure{T}}, H::OE_procedure)    = OE_procedure(H.i1, H.i2,
  convert(T, H.c), convert(T, H.s))

ctranspose(H::OE_procedure) = H
transpose(H::OE_procedure)  = H

function oe_algorithm{T<:AbstractFloat}(x::T, y::T)
  @assert abs(x)  > abs(y) string("oe_Algorithm: |x| > |y| required", abs(x) ," abs(y): " , abs(y))
  if y == zero(T)
    s = zero(T)
    c = one(T)
  else
    s = y/x
    c = sign(x)*sqrt(one(T)-s)*sqrt(one(T)+s)  #sign(x)*
    x = c*x
  end
  return c,s,x
end

"""
    oe_procedure{T}(f::T, g::T, i1::Integer, i2::Integer) -> (H::OE_procedure, r::T)
Computes the Hyperbolic rotation `H` such that for any vector `x` where
```
x[i1] = f
x[i2] = g
```
the result of the multiplication
```
y = H*x
```
has the property that
```
y[i1] = r
y[i2] = 0
```

# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""
function oe_procedure{T}(f::T, g::T, i1::Integer, i2::Integer)
  if i1 == i2
    throw(ArgumentError("Indices must be distinct."))
  end
  if f == g == zero(T)
    throw(ArgumentError("f and g can't both be zero."))
  elseif abs(f) > abs(g) # switch indices
    c,s,x = oe_algorithm(f, g)
    return OE_procedure(i1,i2,convert(T,c),convert(T,s)), x
  else
    c,s,x = oe_algorithm(g, f)
    return OE_procedure(i1,i2,convert(T,c),convert(T,s)), x
  end
end
"""
    oe_procedure(A::AbstractArray, i1::Integer, i2::Integer, j::Integer) -> (H::OE_procedure, r)
Computes the Hyperbolic rotation `H` and scalar `r` such that the result of the multiplication
```
B = H*A
```
has the property that
```
B[i1,j] = r
B[i2,j] = 0
```

# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""
oe_procedure(A::AbstractMatrix, i1::Integer, i2::Integer, j::Integer) =
    oe_procedure(A[i1,j], A[i2,j],i1,i2)


"""
    oe_procedure(A::AbstractVector, i1::Integer, i2::Integer) -> (H::OE_procedure, r)
Computes the Hyperbolic rotation `H` and scalar `r` such that the result of the multiplication
```
B = H*f
```
has the property that
```
B[i1] = r
B[i2] = 0
```
# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""
oe_procedure(A::AbstractVector, i1::Integer, i2::Integer) =
    oe_procedure(A[i1], A[i2], i1, i2)

function getindex(H::OE_procedure, i::Integer, j::Integer)
    if i == j
        if i == H.i1 || i == H.i2
            1/H.c
        else
            one(H.c)
        end
    elseif i == H.i1 && j == H.i2
        -H.s/H.c
    elseif i == H.i2 && j == H.i1
        -H.s/H.c
    else
        zero(H.s)
    end
end

A_mul_B!(H1::OE_procedure, H2::OE_procedure) = error("Operation not supported.")

function A_mul_B!{M1<:AbstractMatrix,T<:AbstractFloat}(H::OE_procedure{T}, A::M1)
  m, n = size(A, 1), size(A, 2)
  if H.i2 > m
    throw(DimensionMismatch("column indices for rotation are outside the matrix"))
  end
  _oe_mul!(view(A,H.i1,:),view(A,H.i2,:),H)
  return A
end

function A_mul_Bc!{M1<:AbstractMatrix,T<:AbstractFloat}(A::M1, H::OE_procedure{T})
    m, n = size(A, 1), size(A, 2)
    if H.i2 > n
        throw(DimensionMismatch("row indices for rotation are outside the matrix"))
    end
    _oe_mul!(view(A,:,H.i1),view(A,:,H.i2),H)
    #_oe_mul!(A,H)
    return A
end

function _oe_mul!{V1<:AbstractVector,T2<:AbstractFloat}(x::V1,y::V1, H::OE_procedure{T2})
  T = promote_type(eltype(x),T2,Float32)
  # abs(x) > abs(y) assumed
  c,s = H.c,H.s
  scale!(x, one(T)/c)
  axpy!(-s/c, y, x)
  scale!(c, y)
  axpy!(-s, x, y)
end
