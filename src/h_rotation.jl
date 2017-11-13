

"""
    H_procedure(i1,i2,α,β,ρ) -> H
A hyperbolic rotation operator implemented by the H-procedure.

The fields `α`, `β` and `ρ = β/α` represent the hyperbolic rotation.

The `H_procedure` type supports left multiplication `H*A` and conjugated transpose
right multiplication `A*H'`. The type doesn't have a `size` and can therefore
be multiplied with matrices of arbitrary size as long as `i2<=size(A,2)` for
`H*A` or `i2<=size(A,1)` for `A*H`.
"""
immutable H_procedure{T}
  i1::Int
  i2::Int
  α::T
  β::T
  ρ::T
  c::Vector{T}

  @compat function (::Type{H_procedure}){T}(i1::Int, i2::Int, α::T, β::T, ρ::T)
    c    = Array{T}(3)
    R    = sqrt(α-β)*sqrt(α+β)    # sqrt(α^2-β^2)
    c[1] = α/R
    c[2] = (α+β)/R
    c[3] = (abs(α)- abs(β))/abs(α)        # 1-abs(ρ)
    new{T}(i1,i2,α,β,ρ,c)
  end

  @compat function (::Type{H_procedure}){T}(i1::Int, i2::Int, α::T, β::T, ρ::T, c::Array{T})
    new{T}(i1,i2,α,β,ρ,c)
  end
end

convert{T}(::Type{H_procedure{T}}, H::H_procedure{T}) = H
convert{T}(::Type{H_procedure{T}}, H::H_procedure)    = H_procedure(H.i1, H.i2,
  convert(T, H.α), convert(T, H.β), convert(T, H.ρ), convert(Vector{T}, H.c))

ctranspose(H::H_procedure) = H
transpose(H::H_procedure)  = H

function h_Algorithm{T}(x::T, y::T)
  @assert abs(x)  > abs(y) string("h_Algorithm: |x| > |y| required", abs(x) ," abs(y): " , abs(y))
  c    = Array{T}(3)
  ρ    = y/x
  tmp  = sqrt(1-ρ)*sqrt(1+ρ)
  α    = abs(x)*tmp
  β    = α*ρ
  R    = sign(x)*sqrt(α-β)*sqrt(α+β)    # this sign(x) is motivated by the sign in
  c[1] = α/R                            # the corresponding Slicot routine which uses oe-procedure
  c[2] = (α+β)/R
  c[3] = (abs(α)- abs(β))/abs(α)
  return α, β, ρ, c
end

"""
    h_procedure{T}(f::T, g::T, i1::Integer, i2::Integer) -> (H::H_procedure, r::T)
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
function h_procedure{T}(f::T, g::T, i1::Integer, i2::Integer)
  if i1 == i2
    throw(ArgumentError("Indices must be distinct."))
  end
  if f == g == zero(T)
    throw(ArgumentError("f and g can't both be zero."))
  elseif abs(f) > abs(g) # switch indices
    α,β,ρ,c = h_Algorithm(f, g)
    return H_procedure(i1,i2,convert(T,α),convert(T,β),convert(T,ρ),convert(Vector{T},c)), α
  else
    α,β,ρ,c = h_Algorithm(g, f)
    return H_procedure(i1,i2,convert(T,α),convert(T,β),convert(T,ρ),convert(Vector{T},c)), α
  end
end
"""
    h_procedure(A::AbstractArray, i1::Integer, i2::Integer, j::Integer) -> (H::H_procedure, r)
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
h_procedure(A::AbstractMatrix, i1::Integer, i2::Integer, j::Integer) =
    h_procedure(A[i1,j], A[i2,j],i1,i2)


"""
    h_procedure(A::AbstractVector, i1::Integer, i2::Integer) -> (H::H_procedure, r)
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
h_procedure(A::AbstractVector, i1::Integer, i2::Integer) =
    h_procedure(A[i1], A[i2], i1, i2)

function getindex(H::H_procedure, i::Integer, j::Integer)
    if i == j
        if i == H.i1 || i == H.i2
            H.c[1]
        else
            one(H.α)
        end
    elseif i == H.i1 && j == H.i2
        -H.c[1]*H.ρ
    elseif i == H.i2 && j == H.i1
        -H.c[1]*H.ρ
    else
        zero(H.α)
    end
end

A_mul_B!(H1::H_procedure, H2::H_procedure) = error("Operation not supported.")

function A_mul_B!{M1<:AbstractMatrix,T<:AbstractFloat}(H::H_procedure{T}, A::M1)
  m, n = size(A, 1), size(A, 2)
  if H.i2 > m
    throw(DimensionMismatch("column indices for rotation are outside the matrix"))
  end
  @inbounds @simd for i = 1:n
    x, y = A[H.i1,i], A[H.i2,i]
    if abs(x) > abs(y)
      x1,y1 = _h_mul(x,y,H)
    else
      y1,x1 = _h_mul(y,x,H)
    end
    A[H.i1,i] = x1
    A[H.i2,i] = y1
  end
  return A
end

function A_mul_Bc!{M1<:AbstractMatrix,T<:AbstractFloat}(A::M1, H::H_procedure{T})
    m, n = size(A, 1), size(A, 2)
    if H.i2 > n
        throw(DimensionMismatch("row indices for rotation are outside the matrix"))
    end
    @inbounds @simd for i = 1:m
        x, y = A[i,H.i1], A[i,H.i2]
        if abs(x) > abs(y)
          x1,y1 = _h_mul(x,y,H)
        else
          y1,x1 = _h_mul(y,x,H)
        end
        A[i,H.i1] = x1
        A[i,H.i2] = y1
    end
    return A
end

function _h_mul{T1<:AbstractFloat,T2<:AbstractFloat}(x::T1, y::T1, H::H_procedure{T2})
  # abs(x) > abs(y) assumed
  α,β,ρ    = H.α,  H.β,  H.ρ
  c        = H.c
  f        = (β*y)/(α*x)
  if f < 0.5
    xi = 1 - f
  else
    d2 = (abs(x)-abs(y))/abs(x)
    xi = c[3] + d2 - c[3]*d2
  end
  x1 = x*xi*c[1]
  y1 = x1 - c[2]*(x-y)
  # make sure that abs(x1) > abs(y1)
  y1 = abs(x1) < abs(y1) ? abs(x1)*(1-3*eps(T1))*sign(y1) :
                           y1
  return x1,y1
end
