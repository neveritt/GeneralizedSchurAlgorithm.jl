println("Starting hyperbolic rotations test...")

# Apply rotations h-procedure and oe_rotation to a random matrix and compare to
# reference implementation
ϵ_t  = 1e-10
M    = 10
N    = 30
C    = randn(N,M)
Ch   = copy(C)
Coe  = copy(C)
Cr   = copy(C)
idx2 = 5
for idx1 = 1:4
  # make sure rotation is defined properly
  if abs(Ch[idx1,idx1]) < abs(Ch[idx2,idx1])
    Ch[idx1,idx1],  Ch[idx2,idx1]  = Ch[idx2,idx1],  Ch[idx1,idx1]
    Coe[idx1,idx1], Coe[idx2,idx1] = Coe[idx2,idx1], Coe[idx1,idx1]
    Cr[idx1,idx1],  Cr[idx2,idx1]  = Cr[idx2,idx1],  Cr[idx1,idx1]
  end
  # h-procedure
  h1 = h_procedure(Ch[idx1,idx1],Ch[idx2,idx1],idx1,idx2)[1]
  A_mul_B!(h1, Ch)
  @test norm(Ch[idx2,idx1]) < ϵ_t

  # oe-procedure
  h2 = oe_procedure(Coe[idx1,idx1],Coe[idx2,idx1],idx1,idx2)[1]
  A_mul_B!(h2, Coe)
  @test norm(Coe[idx2,idx1]) < ϵ_t

  # reference implementation
  ρ = Cr[idx2,idx1]/Cr[idx1,idx1]
  c = 1/sqrt(1-ρ^2)
  s = -ρ*c
  for i = 1:M
    x1, y1 = Cr[idx1,i],Cr[idx2,i]
    Cr[idx1,i] = c*x1 + s*y1
    Cr[idx2,i] = c*y1 + s*x1
  end
end

# oe and h has a sign difference in their rotations
@test norm(abs(Ch)-abs(Cr))  < ϵ_t
@test norm(abs(Coe)-abs(Cr)) < ϵ_t

# right-side multiplication
C    = randn(N,M)
Ch   = copy(C)
Coe  = copy(C)
Cr   = copy(C)
idx2 = 5
for idx1 = 1:4
  # make sure rotation is defined properly
  if abs(Ch[idx1,idx1]) < abs(Ch[idx1,idx2])
    Ch[idx1,idx1],  Ch[idx1,idx2]  = Ch[idx1,idx2],  Ch[idx1,idx1]
    Coe[idx1,idx1], Coe[idx1,idx2] = Coe[idx1,idx2], Coe[idx1,idx1]
    Cr[idx1,idx1],  Cr[idx1,idx2]  = Cr[idx1,idx2],  Cr[idx1,idx1]
  end
  # h-procedure
  h1 = h_procedure(Ch[idx1,idx1],Ch[idx1,idx2],idx1,idx2)[1]
  A_mul_Bc!(Ch, h1)
  @test norm(Ch[idx1,idx2]) < ϵ_t

  # oe-procedure
  h2 = oe_procedure(Coe[idx1,idx1],Coe[idx1,idx2],idx1,idx2)[1]
  A_mul_Bc!(Coe, h2)
  @test norm(Coe[idx1,idx2]) < ϵ_t

  # reference implementation
  ρ = Cr[idx1,idx2]/Cr[idx1,idx1]
  c = 1/sqrt(1-ρ^2)
  s = -ρ*c
  for i = 1:N
    x1, y1 = Cr[i,idx1],Cr[i,idx2]
    Cr[i,idx1] = c*x1 + s*y1
    Cr[i,idx2] = c*y1 + s*x1
  end
end

# oe and h has a sign difference in their rotations
@test norm(abs(Ch)-abs(Cr))  < ϵ_t
@test norm(abs(Coe)-abs(Cr)) < ϵ_t
