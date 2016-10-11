function gen_toeplitz(N=1000, m=200, σ = 1, λ=1e1, n=3)
  # model orders
  na, nb = n, n
  nk = 1

  # Define the true system
  a = [0.5, 0.2, 0.1]
  b = [0.1, 0.3, 0.2]

  # generate input data+noise and simulate output
  A = [1;a]
  B = [0;b]

  u = σ*randn(N)
  e = sqrt(λ)*randn(N)
  y = filt(B,A,u) + filt(1,A,e)

  # generate Toeplitz matrix
  col = hcat(-y[1:end-1],u[1:end-1])
  row = zeros(1,m)
  row[1,1:2] = col[1,1:2]
  Y = y[2:end]
  T = Toeplitz(col,row)

  return T,Y,a,b
end
