      SUBROUTINE RK3_TI_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)
C     ***********************************************************************
C
C       RK3_TI_SOLVE integrates a scalar stochastic differential equation (SDE)
C       using a third-order, time-invariant Runge-Kutta method.
C
C       The SDE has the form:
C
C            dX(t) = F(X) dt + Q * G(X) dW(t)
C
C       Numerical Integration (RK3 step):
C
C            X_{n+1} = X1 + A41*K1 + A42*K2 + A43*K3
C
C            with coefficients:
C              A21 = 1.52880952525675
C              A31 = 0.0, A32 = 0.51578733443615
C              A41 = 0.53289582961739, A42 = 0.25574324768195
C              A43 = 0.21136092270067
C              Q1  = 1.87653936176981, Q2 = 3.91017166264989,
C              Q3  = 4.73124353935667
C
C       This routine advances X from T0 to TN in N steps.
C
C       License:
C            This code is distributed under the MIT license.
C
C       Author:
C            John Burkardt
C
C       Modified:
C            20 June 2010
C
C       Reference:
C            - Kasdin, Jeremy. "Runge-Kutta algorithm for numerical integration
C              of stochastic differential equations." Journal of Guidance, Control,
C              and Dynamics, 18(1), 114-120, 1995.
C            - Kasdin, Jeremy. "Discrete Simulation of Colored Noise and
C              Stochastic Processes." Proceedings of the IEEE, 83(5), 802-827, 1995.
C
C       Parameters:
C
C         Input, external double precision F
C             Deterministic function F(X) in the SDE.
C
C         Input, external double precision G
C             Stochastic function G(X) in the SDE.
C
C         Output, double precision X(0:N)
C             Array of solution values at each time step.
C
C         Input, double precision T0
C             Initial time.
C
C         Input, double precision TN
C             Final time.
C
C         Input, double precision X0
C             Initial state value.
C
C         Input, integer N
C             Number of time steps.
C
C         Input, double precision Q
C             Spectral density of the driving white noise.
C
C         Input, integer SEED
C             Seed for the random number generator.
C
C     ***********************************************************************

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
      DOUBLE PRECISION A21, A31, A32, A41, A42, A43
      DOUBLE PRECISION Q1, Q2, Q3
      DOUBLE PRECISION K1, K2, K3
      DOUBLE PRECISION W1, W2, W3
      DOUBLE PRECISION T1, T2, T3
      DOUBLE PRECISION X1, X2, X3

      DOUBLE PRECISION TEMP

      A21 = 1.52880952525675D+00
      A31 = 0.0D+00
      A32 = 0.51578733443615D+00
      A41 = 0.53289582961739D+00
      A42 = 0.25574324768195D+00
      A43 = 0.21136092270067D+00

      Q1 = 1.87653936176981D+00
      Q2 = 3.91017166264989D+00
      Q3 = 4.73124353935667D+00

      H = (TN - T0) / DBLE(N)

      X(0) = X0
      T = T0

      TEMP = X(0)

      CALL TABLE_HEADER('RK-3 Time-Invariant')
      CALL TABLE_ROW(0, T, TEMP)

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

            TEMP = X(I)

            CALL TABLE_ROW(I, T, TEMP)

      END DO

      CALL TABLE_FOOTER()

      RETURN
      END
