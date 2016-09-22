# generalized schur algorithm


"""
  schuralgorithm(G, F, p, q) -> L, D

Backward stable [version][1] of generalized Schur Algorithm for indefinite matrices.

# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""
function schuralgorithm(Gin,F,p,q)
  r = size(Gin,1)
  N = size(F,1)
  @assert size(F,1) == size(F,2) string("F needs to be square")
  @assert size(Gin,2) == size(F,1) string("G needs to have same number of rows as F")
  @assert p+q       == r string("p+q needs to equal r")

  L = zeros(Float64,N,N)
  D = zeros(Int,N)
  G = copy(Gin)
  for i = 1:n
    if sumabs2(G[1:p,i]) - sumabs2(G[p+1:end,i]) > 0
      _step(view(G,:,i:N), view(F,i:N,i:N), view(L,i,i:N),p,q,true)
      D[i] = 1
    else
      _step(view(G,:,i:N), view(F,i:N,i:N), view(L,i,i:N),p,q,false)
      D[i] = -1
    end
  end
  L, D
end

function _step{M1<:StridedMatrix, M2<:AbstractMatrix}(
    G::M1,F::M2,L,p::Int,q::Int,pos::Bool)
  @assert eltype(G) == eltype(F) == eltype(L)
  T = eltype(G)
  N = size(G,2)

  # zero out first 2:p column of first row
  _householder_left(view(G,1:p,1:N))

  # zero out first p+2:l column of first row
  _householder_left(view(G,p+1:p+q,1:N))

  (idx1,idx2)  = pos ? (1,p+1) : (p+1,1)
  if abs(G[idx1,1]) < abs(G[idx2,1]) # ensure positive/negative definite step
    G[idx1,1] = abs(G[idx2,1])*(1 + 3*eps(T))*sign(G[1,1])
  end
  h = h_procedure(G[idx1,1],G[idx2,1],idx1,idx2)[1]
  A_mul_B!(h, G)
  L[:] = G[idx1,:]
  G[idx1,:] = (L.'*F).'
end
