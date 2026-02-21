      SUBROUTINE RK1_TV_SOLVE(F, G, T0, TN, X0, N, Q, SEED)
C     *****************************************************************
C
C       RK1_TV_SOLVE solves a scalar stochastic differential equation
C       using a first-order time-varying Runge-Kutta (Euler-Maruyama)
C       scheme with the step formula inlined in this routine.
C
C       The equation has the form:
C
C         dX = F(X,t) dt + Q * G(X,t) dW
C
C       Numerical Method (inlined RK1 step):
C
C         X_{n+1} = X_n + H * F(T_n, X_n) + H * G(T_n, X_n) * W_1
C
C         W_1 = N(0,1) * SQRT(Q/H)
C
C       This routine no longer depends on a separate RK step file.
C
C       Licensing:
C
C         This code is distributed under the MIT license.
C
C     *****************************************************************

      IMPLICIT NONE

C     f2py intent(callback) F, G
C     f2py intent(in) T0, TN, X0, N, Q
C     f2py intent(in,out) SEED
C     f2py intent(out) X(0:N)

      INTEGER N, SEED, I
      DOUBLE PRECISION F, G
      EXTERNAL F, G
            DOUBLE PRECISION R8_NORMAL

      DOUBLE PRECISION T0, TN, T
      DOUBLE PRECISION X0, X(0:N)
      DOUBLE PRECISION H, Q
            DOUBLE PRECISION A21, Q1
            DOUBLE PRECISION K1, W1, T1, X1

            A21 = 1.0D+00
            Q1 = 1.0D+00

      H = ( TN - T0 ) / DBLE(N)

      X(0) = X0
      T = T0

            CALL SDE_TABLE_HEADER('RK-1 Time-Variant')
            CALL SDE_TABLE_ROW(0, T, X(0))

      DO I = 1, N

         T = T0 + H * DBLE(I)

                  T1 = T
                  X1 = X(I-1)
                  W1 = R8_NORMAL(SEED) * SQRT(Q1 * Q / H)
                  K1 = H * F(T1, X1) + H * G(T1, X1) * W1
                  X(I) = X1 + A21 * K1

                  CALL SDE_TABLE_ROW(I, T, X(I))

      END DO

            CALL SDE_TABLE_FOOTER()

      RETURN
      END
