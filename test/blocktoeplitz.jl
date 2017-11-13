println("Starting BlockToeplitz test...")

k = 2
l = 5
p = 20
s = 16
M = k*p
N = l*s
col = randn(M,l)
row = randn(k,N)
colc = randn(M,l)+randn(M,l)im
rowc = randn(k,N)+randn(k,N)im

row[1:k,1:l]  = col[1:k,1:l]
rowc[1:k,1:l] = colc[1:k,1:l]
T = Toeplitz(col,row)
Tc = Toeplitz(colc,rowc)

# construction
@test_throws DomainError Toeplitz(col,col[1:k,1:l]+ones(k,l))

Tf = full(T)
Tcf = full(Tc)

# transpose
@test T.' ≈ Tf.'
@test Tc'  ≈ Tcf'
size(T,3)
# dimensions
@test size(T)           == (M,N)
@test size(T,3)         == 1
@test blocksize(T)      == (p,s)
@test blocksize(T,3)    == 1
@test sizeofblock(T)    == (k,l)
@test sizeofblock(T,3)  == 1

# getblock
@test getblock(T,-p+1) ≈ col[M-k+1:M,:]
@test getblock(T, s-1) ≈ row[:,N-l+1:N]
@test_throws BoundsError getblock(T, s)

# getindex
for i in eachindex(T)
  @test T[i]  ≈ Tf[i]
  @test Tc[i] ≈ Tcf[i]
end
@test_throws BoundsError getindex(T,M,N+1)

# getcol and getrow
@test getcol(T) ≈ col
@test getrow(T) ≈ row

# conversion to full matrix
@test convert(Matrix, T) ≈ Tf

# multiplication
C     = randn(N,M)
Cv    = randn(N)
Cl    = randn(N+1,M)
Cvl   = randn(N+1)
Cc    = complex.(C, randn(N,M))
@test T*C  ≈ Tf*C
@test_throws DimensionMismatch T*Cl
@test T*Cv ≈ Tf*Cv
@test_throws DimensionMismatch T*Cvl
@test Tc*Cc ≈ Tcf*Cc
# works in 0.5
#@test C*T  ≈ C*Tf
#@test Cc*Tc ≈ Cc*Tcf


# At_mul_B
C     = randn(M,N)
Cv    = randn(M)
Cl    = randn(M+1,N)
Cvl   = randn(M+1)
Cc    = complex.(C, randn(M,N))
@test T.'*C   ≈ Tf.'*C
@test_throws DimensionMismatch T.'*Cl
@test T.'*Cv ≈ Tf.'*Cv
@test_throws DimensionMismatch T.'*Cvl
@test Tc.'*Cc ≈ Tcf.'*Cc
# works in 0.5
#@test C*T.'   ≈ C*Tf.'
#@test Cc*Tc.' ≈ Cc*Tcf.'

@test T'*C   ≈ Tf'*C
@test Tc'*Cc ≈ Tcf'*Cc
# works in 0.5
#@test C*T'   ≈ C*Tf'
#@test Cc*Tc' ≈ Cc*Tcf'
