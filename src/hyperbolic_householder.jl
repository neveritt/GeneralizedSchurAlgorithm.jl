function _hyper_householder{M1<:AbstractMatrix}(A::M1,g::Int,p::Int,q::Int,work::M1)
  gp   = g[1:p]
  gq   = g[p+1:p+q]
  gJgT = sum(abs2, gp) - sum(abs2, gq)
  β    = 2/gJgT
  B    = β*[gp*gp' -gp*gq'; gq*gp' -gq*gq']

  T = eltype(A)
  H = eye(T,length(g)) - B

  # perform multiplication in place
  A_mul_B!(work,H,A)
end

function _hyper_householder{M1<:AbstractMatrix}(A::M1,g,p::Int,q::Int)
  gp   = g[1:p]
  gq   = g[p+1:p+q]
  gJgT = sum(abs2, gp) - sum(abs2, gq)
  β    = 2/gJgT
  B    = β*[gp*gp' -gp*gq'; gq*gp' -gq*gq']

  T = eltype(A)
  H = eye(T,length(g)) - B

  A[:,:] = H * A[:,:]
end

# Fast Reliable Algorithms for Matrices with Structure
# T. Kailath and A. H. Sayed
# 1999

function _householder{M1<:AbstractMatrix, M3<:AbstractMatrix}(
  A::M1,j::Int,work::M1,H::M3,w::M3)
  @assert eltype(A) == eltype(work) == eltype(H) == eltype(w)
  @assert size(work) == size(A)            string("work must have correct size")
  @assert size(w) == (size(A,1),1)         string("w must have correct size")
  @assert size(H) == (size(A,1),size(A,1)) string("H must have correct size")
  T = eltype(A)
  w[:,1] = A[:,1]
  σ = norm(w) # σ
  if σ == zero(T)
    return
  else
    β     = 1/(σ*(σ + abs(A[j,1])))
    s     = A[j,1] > 0 ? one(T) : -one(T) # s
    w[j]  = s*(σ + abs(A[j,1]))
    n = length(w)
    H[:,:] = eye(T,n) - β*w*w.'
    A_mul_B!(work,H,A)
    return
  end
end
