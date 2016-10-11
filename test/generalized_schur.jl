println("Starting generalized version test...")

A = gen_toeplitz(1000, 20)[1]
T = eltype(A)

# construct generator G used in qrtoeplitz
m,n = blocksize(A)
k,l = sizeofblock(A)
M,N = size(A)
C,R = qr(getcol(A))
S = A.'*C

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

# construct shift matrix Z
Z1 = spdiagm(ones(T,l*(n-1)),l,l*n,l*n)
Z2 = spdiagm(ones(T,k*(m-1)),k,k*m,k*m)
Z  = [Z1 spzeros(T,n*l,m*k); spzeros(T,m*k,n*l) Z2]
p = k+l
q = k+l
G = G.'

L,D = schuralgorithm(G,Z,p,q,n*l)
Q,R = L[1:n*l, n*l+(1:m*k)].', triu(L[1:n*l, 1:n*l])

@test norm(Q*R-A) < 1e-10

# compare with specialized implementation
Qt,Rt = qrtoeplitz(A)
@test norm(Q*R-Qt*Rt) < 1e-10
