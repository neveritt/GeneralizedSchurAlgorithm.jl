println("Starting Toeplitz test...")

# Define the true system
a = [0.5, 0.2, 0.1]
b = [0.1, 0.3, 0.2]

# model orders
na, nb = 3, 3
nk = 1
n = [na, nb, nk]

# generate input data+noise and simulate output
A = [1;a]
B = [0;b]
N = 3000

u = randn(N)
λ = 1e1
e = sqrt(λ)*randn(N)
y = filt(B,A,u) + filt(1,A,e)

col = hcat(-y[1:end-1],u[1:end-1])
row = zeros(1,6)
row[1,1:2] = col[1,1:2]
Y = y[2:end]
T = Toeplitz(col,row)

# test QR - factorization
Q,R = qrtoeplitz(T)

@test norm(Q*R-T) < 1e-10

# test least squares solver
x̂ = lstoeplitz(T,Y)

Q,R = qr(full(T))
x = R\(Q.'*Y)
@test norm(x̂-x) < 1e-3

# test against true parameters
x0 = zeros(length(a)+length(b))
x0[1:2:end-1] = a
x0[2:2:end]   = b
@test norm(x̂-x0) < 2λ
