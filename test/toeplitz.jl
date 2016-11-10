println("Starting Toeplitz test...")

ϵ_t  = 1e-10

N = 1000; m = 200; σ = 1; λ=1e1
T,Y,a,b = gen_toeplitz(N, m, σ, λ)
na = length(a); nb = length(b)

# test QR - factorization
Q,R = qrtoeplitz(T)
Q2,R2 = qrtoeplitz(getcol(T),getrow(T))
@test norm(Q   - Q2) < ϵ_t
@test norm(R   - R2) < ϵ_t
@test norm(Q*R -  T) < ϵ_t

# test least squares solver
x̂  = lstoeplitz(T,Y)
x1 = lstoeplitz(getcol(T),getrow(T),Y)
@test norm(x̂-x1) < ϵ_t

Q,R = qr(full(T))
x = R\(Q.'*Y)
@test norm(x̂-x) < 1e-3

# test against true parameters
x0 = zeros(size(T,2)) #zeros(length(a)+length(b))
x0[1:2:na+nb-1] = a
x0[2:2:na+nb]   = b
@test norm(x̂-x0) < 5λ/σ^2*m/N
