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

Tf = full(T)
Tcf = full(Tc)

# transpose
@test T.' ≈ Tf.'
@test Tc'  ≈ Tcf'

# dimensions
@test size(T)        == (M,N)
@test blocksize(T)   == (p,s)
@test sizeofblock(T) == (k,l)

# getblock
@test getblock(T,-p+1) ≈ col[M-k+1:M,:]
@test getblock(T, s-1) ≈ row[:,N-l+1:N]

# getindex
for i in eachindex(T)
  @test T[i]  ≈ Tf[i]
  @test Tc[i] ≈ Tcf[i]
end

# multiplication
C  = randn(N,M)
Cc = complex(C, randn(N,M))
@test T*C  ≈ Tf*C
@test Tc*Cc ≈ Tcf*Cc
# works in 0.5
#@test C*T  ≈ C*Tf
#@test Cc*Tc ≈ Cc*Tcf

# At_mul_B
C = randn(M,N)
Cc = complex(C, randn(M,N))
@test T.'*C   ≈ Tf.'*C
@test Tc.'*Cc ≈ Tcf.'*Cc
# works in 0.5
#@test C*T.'   ≈ C*Tf.'
#@test Cc*Tc.' ≈ Cc*Tcf.'

@test T'*C   ≈ Tf'*C
@test Tc'*Cc ≈ Tcf'*Cc
# works in 0.5
#@test C*T'   ≈ C*Tf'
#@test Cc*Tc' ≈ Cc*Tcf'
