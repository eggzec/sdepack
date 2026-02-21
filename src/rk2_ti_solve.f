      SUBROUTINE RK2_TI_SOLVE(F, G, T0, TN, X0, N, Q, SEED)
C     *****************************************************************
C
C       RK2_TI_SOLVE solves a scalar stochastic differential equation
C       using a second-order time-invariant Runge-Kutta method with
C       the step formula inlined in this routine.
C
C       The equation has the form:
C
C         dX = F(X) dt + Q * G(X) dW
C
C       Numerical Method (inlined RK2 step):
C
C         K1 = H * F(X1) + H * G(X1) * W1
C         K2 = H * F(X2) + H * G(X2) * W2
C
C         X2     = X1 + A21 * K1
C         X_{n+1}= X1 + A31 * K1 + A32 * K2
C
C         A21 = 1, A31 = 0.5, A32 = 0.5, Q1 = Q2 = 2
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
            DOUBLE PRECISION A21, A31, A32
            DOUBLE PRECISION Q1, Q2
            DOUBLE PRECISION K1, K2, W1, W2
            DOUBLE PRECISION T1, T2
            DOUBLE PRECISION X1, X2

            A21 = 1.0D+00
            A31 = 0.5D+00
            A32 = 0.5D+00

            Q1 = 2.0D+00
            Q2 = 2.0D+00

      H = ( TN - T0 ) / DBLE(N)

      X(0) = X0
      T = T0

            CALL SDE_TABLE_HEADER('RK-2 Time-Invariant')
            CALL SDE_TABLE_ROW(0, T, X(0))

      DO I = 1, N

         T = T0 + H * DBLE(I)

                  T1 = T
                  X1 = X(I-1)
                  W1 = R8_NORMAL(SEED) * SQRT(Q1 * Q / H)
                  K1 = H * F(X1) + H * G(X1) * W1

                  T2 = T1 + A21 * H
                  X2 = X1 + A21 * K1
                  W2 = R8_NORMAL(SEED) * SQRT(Q2 * Q / H)
                  K2 = H * F(X2) + H * G(X2) * W2

                  X(I) = X1 + A31 * K1 + A32 * K2

                  CALL SDE_TABLE_ROW(I, T, X(I))

      END DO

            CALL SDE_TABLE_FOOTER()

      RETURN
      END
