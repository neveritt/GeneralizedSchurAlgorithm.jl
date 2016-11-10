println("Starting generalized version test...")

ϵ_t  = 1e-10

A,Y = gen_toeplitz(1000, 20)
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

@test norm(Q*R-A) < ϵ_t

# compare with specialized implementation
Qt,Rt = qrtoeplitz(A)
@test norm(Q*R-Qt*Rt) < ϵ_t

# construct generator G used in lstoeplitz
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

β = 4*(m*k*n*l)^(4/3)*eps(T) # tuning parameter to ensure negative steps
GQ[:,1:l]       = C
GQ[:,l+k+(1:l)] = C
for i = 1:k
  GQ[i,l+i] = one(T)
  GQ[i,2k+2l+i] = 1+β
end
Z1 = spdiagm(ones(T,l*(n-1)),l,l*n,l*n)
Z2 = spdiagm(ones(T,k*(m-1)),k,k*m,k*m)
Z  = [Z1 spzeros(T,n*l,m*k); spzeros(T,m*k,n*l) Z2]

# iterate
p = k+l
q = k+l+k
G = G.'
N = m*k+n*l
L,D   = schuralgorithm(G,Z,p,q,N)

Q, R  = L[1:n*l, n*l+(1:m*k)].', triu(L[1:n*l, 1:n*l])
D     = tril(L[n*l+(1:m*k), n*l+(1:m*k)].')
x     = R\(D\Q).'*(D\Y)

# compare with specialized implementation
xt = lstoeplitz(A,Y)
@test norm(x-xt) < ϵ_t
