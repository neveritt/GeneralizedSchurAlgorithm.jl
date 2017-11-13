# generalized schur algorithm

"""
  schuralgorithm(G, F, p, q) -> L, D

Backward stable [version][1] of generalized Schur Algorithm for indefinite
matrices. The algorithm modifies the generator G.

# References

-  [1]: T. Kailath and A. H. Sayed, Fast Reliable Algorithms for Matrices
      with Structure, Society for Industrial and Applied Mathematics, 1999.
"""
function schuralgorithm(G,F,p,q,n)
  r = size(G,1)
  N = size(F,1)
  @assert size(F,1) == size(F,2) "F needs to be square"
  @assert size(G,2) == size(F,1) "G needs to have same number of rows as F"
  @assert p+q       == r         "p+q needs to equal r"

  L = zeros(Float64,n,N)
  D = zeros(Int,n)
  for i = 1:n
    if sum(abs2, G[1:p,i]) - sum(abs2, G[p+1:end,i]) > 0
      _step(view(G,:,i:N), view(L,i,i:N),p,q,true)
      D[i] = 1
      G[1,i:N] = (L[i:i,i:N]*F[i:N,i:N]).'
    else
      _step(view(G,:,i:N), view(L,i,i:N),p,q,false)
      D[i] = -1
      G[p+1,i:N] = (L[i:i,i:N]*F[i:N,i:N]).'
    end
  end
  L, D
end

function _step{M1<:StridedMatrix}(
    G::M1,L,p::Int,q::Int,pos::Bool)
  @assert eltype(G) == eltype(L)
  T = eltype(G)
  N = size(G,2)

  # zero out row 2:p of first column
  _householder_left(view(G,1:p,1:N))

  # zero out row p+2:p+q of first column
  _householder_left(view(G,p+1:p+q,1:N))

  (idx1,idx2)  = pos ? (1,p+1) : (p+1,1)
  if abs(G[idx1,1]) < abs(G[idx2,1]) # ensure positive/negative definite step
    G[idx1,1] = abs(G[idx2,1])*(1 + 3*eps(T))*sign(G[1,1])
  end
  h = h_procedure(G[idx1,1],G[idx2,1],idx1,idx2)[1]
  A_mul_B!(h, G)
  L[:] = G[idx1,:]
end
