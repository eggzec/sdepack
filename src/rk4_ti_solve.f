      SUBROUTINE RK4_TI_SOLVE(F, G, T0, TN, X0, N, Q, SEED)
C     *****************************************************************
C
C       RK4_TI_SOLVE solves a scalar stochastic differential equation
C       using a fourth-order time-invariant Runge-Kutta method with
C       the step formula inlined in this routine.
C
C       The equation has the form:
C
C         dX = F(X) dt + Q * G(X) dW
C
C       Numerical Method (inlined RK4 step):
C
C         Four-stage stochastic RK update
C
C           X_{n+1} = X1 + A51*K1 + A52*K2 + A53*K3 + A54*K4
C
C         using Kasdin coefficients for time-invariant systems.
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
            DOUBLE PRECISION A21, A31, A32, A41, A42, A43
            DOUBLE PRECISION A51, A52, A53, A54
            DOUBLE PRECISION Q1, Q2, Q3, Q4
            DOUBLE PRECISION K1, K2, K3, K4
            DOUBLE PRECISION W1, W2, W3, W4
            DOUBLE PRECISION T1, T2, T3, T4
            DOUBLE PRECISION X1, X2, X3, X4

            A21 = 2.71644396264860D+00
            A31 = -6.95653259006152D+00
            A32 = 0.78313689457981D+00
            A41 = 0.0D+00
            A42 = 0.48257353309214D+00
            A43 = 0.26171080165848D+00
            A51 = 0.47012396888046D+00
            A52 = 0.36597075368373D+00
            A53 = 0.08906615686702D+00
            A54 = 0.07483912056879D+00

            Q1 = 2.12709852335625D+00
            Q2 = 2.73245878238737D+00
            Q3 = 11.22760917474960D+00
            Q4 = 13.36199560336697D+00

      H = ( TN - T0 ) / DBLE(N)

      X(0) = X0
      T = T0

            CALL SDE_TABLE_HEADER('RK-4 Time-Invariant')
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

                  T3 = T1 + A31 * H + A32 * H
                  X3 = X1 + A31 * K1 + A32 * K2
                  W3 = R8_NORMAL(SEED) * SQRT(Q3 * Q / H)
                  K3 = H * F(X3) + H * G(X3) * W3

                  T4 = T1 + A41 * H + A42 * H + A43 * H
                  X4 = X1 + A41 * K1 + A42 * K2 + A43 * K3
                  W4 = R8_NORMAL(SEED) * SQRT(Q4 * Q / H)
                  K4 = H * F(X4) + H * G(X4) * W4

                  X(I) = X1 + A51 * K1 + A52 * K2 + A53 * K3 + A54 * K4

                  CALL SDE_TABLE_ROW(I, T, X(I))

      END DO

            CALL SDE_TABLE_FOOTER()

      RETURN
      END
