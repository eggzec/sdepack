      SUBROUTINE RK1_TI_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)
C     *****************************************************************
C
C       RK1_TI_SOLVE solves a scalar stochastic differential equation
C       using a first-order time-invariant Runge-Kutta (Euler-Maruyama)
C       scheme with the step formula inlined in this routine.
C
C       The equation has the form:
C
C         dX = F(X) dt + Q * G(X) dW
C
C       Numerical Method (inlined RK1 step):
C
C         X_{n+1} = X_n + H * F(X_n) + H * G(X_n) * W_1
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

      INTEGER N
      INTEGER SEED
      INTEGER I

      EXTERNAL F
      EXTERNAL G
      DOUBLE PRECISION F
      DOUBLE PRECISION G
      DOUBLE PRECISION R8_NORMAL

      DOUBLE PRECISION T0, TN, T
      DOUBLE PRECISION X0
      DOUBLE PRECISION X(0:N)
      DOUBLE PRECISION H
      DOUBLE PRECISION Q
      DOUBLE PRECISION A21, Q1
      DOUBLE PRECISION K1, W1, X1

      DOUBLE PRECISION TEMP

      A21 = 1.0D+00
      Q1 = 1.0D+00

      H = (TN - T0) / DBLE(N)

      X(0) = X0
      T = T0

      TEMP = X(0)

      CALL TABLE_HEADER('RK-1 Time-Invariant')
      CALL TABLE_ROW(0, T, TEMP)

      DO I = 1, N
            T = T0 + H * DBLE(I)

            X1 = X(I-1)
            W1 = R8_NORMAL(SEED) * SQRT(Q1 * Q / H)
            K1 = H * F(X1) + H * G(X1) * W1
            X(I) = X1 + A21 * K1

            TEMP = X(I)

            CALL TABLE_ROW(I, T, TEMP)

      END DO

      CALL TABLE_FOOTER()

      RETURN
      END
