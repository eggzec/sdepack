      SUBROUTINE RK3_TI_SOLVE(F, G, T0, TN, X0, N, Q, SEED)
C     *****************************************************************
C
C       RK3_TI_SOLVE solves a scalar stochastic differential equation
C       using a third-order time-invariant Runge-Kutta method with
C       the step formula inlined in this routine.
C
C       The equation has the form:
C
C         dX = F(X) dt + Q * G(X) dW
C
C       Numerical Method (inlined RK3 step):
C
C         Three-stage stochastic RK update
C
C           X_{n+1} = X1 + A41*K1 + A42*K2 + A43*K3
C
C         with coefficients:
C
C           A21 = 1.52880952525675
C           A31 = 0.0, A32 = 0.51578733443615
C           A41 = 0.53289582961739, A42 = 0.25574324768195
C           A43 = 0.21136092270067
C           Q1  = 1.87653936176981, Q2 = 3.91017166264989
C           Q3  = 4.73124353935667
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
            DOUBLE PRECISION Q1, Q2, Q3
            DOUBLE PRECISION K1, K2, K3
            DOUBLE PRECISION W1, W2, W3
            DOUBLE PRECISION T1, T2, T3
            DOUBLE PRECISION X1, X2, X3

            A21 = 1.52880952525675D+00
            A31 = 0.0D+00
            A32 = 0.51578733443615D+00
            A41 = 0.53289582961739D+00
            A42 = 0.25574324768195D+00
            A43 = 0.21136092270067D+00

            Q1 = 1.87653936176981D+00
            Q2 = 3.91017166264989D+00
            Q3 = 4.73124353935667D+00

      H = ( TN - T0 ) / DBLE(N)

      X(0) = X0
      T = T0

            CALL SDE_TABLE_HEADER('RK-3 Time-Invariant')
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

                  X(I) = X1 + A41 * K1 + A42 * K2 + A43 * K3

                  CALL SDE_TABLE_ROW(I, T, X(I))

      END DO

            CALL SDE_TABLE_FOOTER()

      RETURN
      END
