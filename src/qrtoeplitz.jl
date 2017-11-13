"""
  qrtoeplitz(col, row) -> Q, R

Fast Q,R factorization of a Toeplitz matrix. The backward stable [version][1] of
the generalized Schur Algorithm is used.

# Examples
```julia
julia> col = reshape([1.0, 2.4, 5.0, 6.1],4,1);

julia> row = [1.0 -5.1 1.2];

julia> Q,R = qrtoeplitz(col,row);
```

# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""

function qrtoeplitz(vc, vr)
  qrtoeplitz(Toeplitz(vc,vr))
end

#
# col = C*R  - qr factorization of first column of A
# S   = A.'*C
# GR =
#  [S_0         0       0          0
#  S_1       A_{1}     S1      A_{-m+1}
#   ⋮          ⋮       ⋮        ⋮
#  S_{n-1}  A_{n-1}  S_{n-1}  A_{-m+n-1}]
#
# GQ =
#  [C_0      I_k     C0     0
#   C_1       0      C1     0
#   ⋮        ⋮      ⋮     ⋮
#   C_{m-1}   0    C_{m-1}  0]
#
function qrtoeplitz{T<:Number}(A::BlockToeplitz{T})
  m,n = blocksize(A)
  k,l = sizeofblock(A)
  M,N = size(A)
  C,R = qr(getcol(A))
  S   = A.'*C

  # construct generator G
  G  = zeros(T,N+M, l+k+l+k)
  GR = view(G, 1:N, 1:2l+2k)
  GQ = view(G, N+1:N+M, 1:2l+2k)
  GR[:,1:l]         = S
  GR[:,k+l+(1:l)]   = S
  GR[1:l,k+l+(1:l)] = zeros(T,l,l)
  for i = 1:n-1
    GR[i*l+(1:l),l+(1:k)]     = getblock(A,i).'
    GR[i*l+(1:l),l+k+l+(1:k)] = getblock(A,i-m).'
  end
  GQ[:,1:l]       = C
  GQ[:,l+k+(1:l)] = C
  for i = 1:k
    GQ[i,l+i] = one(T)
  end

  # iterate
  p = k+l
  q = k+l
  L = Array{T}(n*l,m*k+n*l)
  G = G.'
  N = size(G,2)
  for i = 1:n*l
    @inbounds _step(view(G,:,i:N), view(L,i,i:N),p,q,true)
    _downshift!(view(G,1,:),view(L,i,:),l,k,i,N,l*n)
  end

  return L[1:n*l, n*l+(1:m*k)].', triu(L[1:n*l, 1:n*l])
end
