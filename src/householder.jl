function _householder_left{M1<:StridedMatrix}(A::M1,first::Bool=true)
  T = eltype(A)
  N = size(A,2)
  x = A[:,1]
  n = length(x)
  α = x[1,1]

  β, τ   = larfg!(α, view(x,2:n))
  x[1,1] = one(T)

  A[2:n,1] = zeros(T,n-1)
  A[1,1]   = β[]

  larfx!('L', view(A,:,2:N), x, τ[])
end

for (larf, larfg, larfx, elty) in
  ((:dlarf_,:dlarfg_,:dlarfx_,:Float64),
   (:slarf_,:slarfg_,:slarfx_,:Float32),
   (:zlarf_,:zlarfg_,:zlarfx_,:Complex128),
   (:clarf_,:clarfg_,:clarfx_,:Complex64))

  @eval begin
    # SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
    # *   .. Scalar Arguments ..
    #     CHARACTER          SIDE
    #     INTEGER            LDC, M, N
    #     DOUBLE PRECISION   TAU
    # *   .. Array Arguments ..
    #     DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
    function larf!(side::Char, C::StridedMatrix{$elty},
                    v::StridedVector{$elty}, τ::$elty)
      m,n = size(C)

      ldc = stride(C,2)
      incv = 1
      if side == 'L'
        work = Array{$elty}(n)
      else # side == 'R'
        work = Array{$elty}(m)
      end
      ccall((Base.LinAlg.BLAS.@blasfunc($larf),Base.liblapack_name), Void,
            (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
             Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}),
            &side, &m, &n, v, &incv, &τ, C, &ldc, work)
    end

    # SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
    # *   .. Scalar Arguments ..
    #      INTEGER            INCX, N
    #      DOUBLE PRECISION   ALPHA, TAU
    # *   .. Array Arguments ..
    #     DOUBLE PRECISION   X( * )
    function larfg!(α::$elty, x::StridedVector{$elty})
      n    = length(x)+1
      incx = one(Int)
      τ    = Ref{$elty}()
      β    = Ref{$elty}(α)
      ccall((Base.LinAlg.BLAS.@blasfunc($larfg),Base.liblapack_name), Void,
            (Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}),
            &n, β, x, &incx, τ)

      β, τ
    end

    # SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
    # *   .. Scalar Arguments ..
    #     CHARACTER          SIDE
    #     INTEGER            LDC, M, N
    #     DOUBLE PRECISION   TAU
    # *   .. Array Arguments ..
    #     DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
    function larfx!(side::Char, C::StridedMatrix{$elty},
                    v::StridedVector{$elty}, τ::$elty)
      m,n = size(C)

      ldc = stride(C,2)
      if side == 'L' && m > 10
        work = Array{$elty}(n)
      elseif m > 10 # side == 'R'
        work = Array{$elty}(m)
      else # work is not referenced if H has order < 11
        work = Array{$elty}(0)
      end
      ccall((Base.LinAlg.BLAS.@blasfunc($larfx),Base.liblapack_name), Void,
            (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty},
             Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}),
            &side, &m, &n, v, &τ, C, &ldc, work)
    end
  end
end
