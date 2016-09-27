function lstoeplitz(col, row, b)
  lstoeplitz(BlockToeplitz(col, row), b)
end


# G = [GR; GQ]
#
# GR =
#  [S_0         0       0          0       0
#  S_1       A_{1}     S1      A_{-m+1}    0
#   ⋮          ⋮       ⋮        ⋮       ⋮
#  S_{n-1}  A_{n-1}  S_{n-1}  A_{-m+n-1}]  0]
#
# GR =
#  [C_0      I_k     C0     0   I_k
#   C_1       0      C1     0    0
#   ⋮        ⋮      ⋮     ⋮   ⋮
#   C_{m-1}   0    C_{m-1}  0    0 ]
#
# col = C*R  - qr factorization of first column of A
# S   = A.'*C
"""
  lstoeplitz(col, row, b) -> x, Q, R, D

Backward stable solution of least squares problem Tx = b, where T is a Toeplitz
matrix defined by the first block column `col` and first block row 'row'.

If x solves Tx = b, then it also satisfies

M = [T^TT T^T;
      T    0 ]

M [x; -b] = [0; b]

The backward stable [version][1] of the generalized Schur Algorithm is used to
compute the cholesky factor `L` of `M`

L̂ = [R̂^T  0;
     Q̂    D]

and then computes
x̂ = R̂\(D\Q̂).'(D\b)

# Examples
```julia
julia> col = reshape([1.0, 2.4, 5.0, 6.1],4,1);

julia> row = [1.0 -5.1 1.2];

julia> T   = BlockToeplitz(col,row);

julia> x   = [3.1,2.3,4.2]

julia> b = T*x + 0.1*randn(size(T,1));

julia> x̂,Q,R,D = lstoeplitz(col,row,b);

julia> norm(x-x̂)
```

# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""
function lstoeplitz{T1<:AbstractFloat, M1<:AbstractArray}(A::BlockToeplitz{T1}, b::M1)
  T   = promote_type(T1, eltype(b))
  m,n = blocksize(A)
  k,l = sizeofblock(A)
  M,N = size(A)
  C,R = qr(getcol(A))
  S = A'*C

  # construct generator G
  G = zeros(T,N+M, 2l+3k)
  GR = view(G, 1:N, 1:2l+3k)
  GQ = view(G, N+1:N+M, 1:2l+3k)
  GR[:,1:l]         = S
  GR[:,k+l+(1:l)]   = S
  GR[1:l,k+l+(1:l)] = zeros(T, l,l)
  for i = 1:n-1
    GR[i*l+(1:l),l+(1:k)]     = getblock(A,i).'
    GR[i*l+(1:l),l+k+l+(1:k)] = getblock(A,i-m).'
  end

  β = 4*(m*k*n*l)^(4/3)*eps(Float64)
  GQ[:,1:l]       = C
  GQ[:,l+k+(1:l)] = C
  for i = 1:k
    GQ[i,l+i] = one(T)
    GQ[i,2k+2l+i] = 1+β
  end

  # construct shift matrix Z
  Z1 = spdiagm(ones(l*(n-1)),l,l*n,l*n)
  Z2 = spdiagm(ones(k*(m-1)),k,k*m,k*m)
  Z::SparseMatrixCSC{T,Int} = [Z1 spzeros(T,n*l,m*k); spzeros(T,m*k,n*l) Z2]
  p = k+l
  q = k+l+k
  G = G.'

  L = zeros(T,m*k+n*l,m*k+n*l)
  N = m*k+n*l
  # positive steps
  for i = 1:n*l
    _step(view(G,:,i:N), Z[i:N,i:N], view(L,i,i:N),p,q,true)
  end
  # negative steps
  for i = n*l+1:N
    _step(view(G,:,i:N), Z[i:N,i:N], view(L,i,i:N),p,q,false)
  end

  Q, R  = L[1:n*l, n*l+(1:m*k)].', triu(L[1:n*l, 1:n*l])
  D     = tril(L[n*l+(1:m*k), n*l+(1:m*k)].')
  x     = R\(D\Q).'*(D\b)
  return x, Q, R, D
end
