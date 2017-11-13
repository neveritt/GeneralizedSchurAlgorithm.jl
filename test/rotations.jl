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

# constructors
h1  = h_procedure(C[1,1],Ch[2,1],1,2)[1]
h11 = H_procedure(h1.i1, h1.i2, h1.α, h1.β, h1.ρ)
h12 = h_procedure(C,1,2,1)[1]

for field in fieldnames(h1)
  # H_procedure has a sign difference which is introduce in matrix constructor for increas numerical accuracy
  @test abs.(getfield(h1,field)) ≈ abs.(getfield(h11,field))
  @test getfield(h1,field)      ≈ getfield(h12,field)
end
h2  = oe_procedure(C[1,1],Ch[2,1],1,2)[1]
h21 = OE_procedure(h2.i1, h2.i2, h2.c, h2.s)
h22 = oe_procedure(C,1,2,1)[1]
for field in fieldnames(h2)
  @test getfield(h2,field) ≈ getfield(h21,field) ≈ getfield(h22,field)
end

# test constructor checks
@test_throws ArgumentError h_procedure(C,1,1,1)
@test_throws ArgumentError oe_procedure(C,1,1,1)
@test_throws ArgumentError h_procedure(0,0,1,2)
@test_throws ArgumentError oe_procedure(0,0,1,2)

# test conjugate transpose
@test h1' == h1
@test h2' == h2
@test h1.' == h1
@test h2.' == h2

# left-side multiplication
idx2 = 5
for idx1 = 1:4
  # make sure rotation is defined properly
  if abs.(Ch[idx1,idx1]) < abs.(Ch[idx2,idx1])
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
  c = sign(Cr[idx1,idx1]) * 1/sqrt(1-ρ^2)
  s = -ρ*c
  for i = 1:M
    x1, y1     = Cr[idx1,i],Cr[idx2,i]
    Cr[idx1,i] = c*x1 + s*y1
    Cr[idx2,i] = c*y1 + s*x1
  end

  @test h1[idx1,idx1]   ≈ h2[idx1,idx1] ≈ c
  @test h1[idx2,idx2]   ≈ h2[idx2,idx2] ≈ c
  @test h1[idx1,idx2]   ≈ h2[idx1,idx2] ≈ s
  @test h1[idx2,idx1]   ≈ h2[idx2,idx1] ≈ s
  @test h1[10,10]       ≈ h2[10,10]     ≈ 1
  @test h1[10,1]        ≈ h2[10,1]      ≈ 0
end

@test norm(Ch  - Cr)    < ϵ_t
@test norm(Coe - Cr)    < ϵ_t

# right-side multiplication
C    = randn(N,M)
Ch   = copy(C)
Coe  = copy(C)
Cr   = copy(C)
idx2 = 5
for idx1 = 1:4
  # make sure rotation is defined properly
  if abs.(Ch[idx1,idx1]) < abs.(Ch[idx1,idx2])
    Ch[idx1,idx1],  Ch[idx1,idx2]  = Ch[idx1,idx2],  Ch[idx1,idx1]
    Coe[idx1,idx1], Coe[idx1,idx2] = Coe[idx1,idx2], Coe[idx1,idx1]
    Cr[idx1,idx1],  Cr[idx1,idx2]  = Cr[idx1,idx2],  Cr[idx1,idx1]
  end
  # h-procedure
  h1 = h_procedure(Ch[idx1,idx1],Ch[idx1,idx2],idx1,idx2)[1]
  A_mul_Bc!(Ch, h1)
  @test norm(Ch[idx1,idx2])  < ϵ_t

  # oe-procedure
  h2 = oe_procedure(Coe[idx1,idx1],Coe[idx1,idx2],idx1,idx2)[1]
  A_mul_Bc!(Coe, h2)
  @test norm(Coe[idx1,idx2]) < ϵ_t

  # reference implementation
  ρ = Cr[idx1,idx2]/Cr[idx1,idx1]
  c = sign(Cr[idx1,idx1]) * 1/sqrt(1-ρ^2)
  s = -ρ*c
  for i = 1:N
    x1, y1     = Cr[i,idx1],Cr[i,idx2]
    Cr[i,idx1] = c*x1 + s*y1
    Cr[i,idx2] = c*y1 + s*x1
  end

  @test h1[idx1,idx1]   ≈ h2[idx1,idx1] ≈ c
  @test h1[idx2,idx2]   ≈ h2[idx2,idx2] ≈ c
  @test h1[idx1,idx2]   ≈ h2[idx1,idx2] ≈ s
  @test h1[idx2,idx1]   ≈ h2[idx2,idx1] ≈ s
  @test h1[10,10]       ≈ h2[10,10]     ≈ 1
  @test h1[10,1]        ≈ h2[10,1]      ≈ 0
end

@test norm(Ch  - Cr)    < ϵ_t
@test norm(Coe - Cr)    < ϵ_t
