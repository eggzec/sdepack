      SUBROUTINE RK1_TV_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)
C     ***********************************************************************
C
C       RK1_TV_SOLVE integrates a scalar stochastic differential equation (SDE)
C       using a first-order, time-varying Runge-Kutta method (Euler-Maruyama).
C
C       The SDE has the form:
C
C            dX(t) = F(X,t) dt + Q * G(X,t) dW(t)
C
C       where:
C            - X(t) is the state variable
C            - F(X,t) is the deterministic drift function
C            - G(X,t) is the stochastic diffusion function
C            - Q is the spectral density of the driving white noise
C            - dW(t) is a Wiener process increment
C
C       Numerical Integration (Euler-Maruyama / RK1 step):
C
C            X_{n+1} = X_n + H * F(T_n, X_n) + sqrt(H*Q) * G(T_n, X_n) * Z_n
C
C            where:
C              - H = (TN - T0)/N is the time step
C              - Z_n ~ N(0,1) is a standard normal random variable
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
C             Deterministic function F(X,T) in the SDE.
C
C         Input, external double precision G
C             Stochastic function G(X,T) in the SDE.
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
      DOUBLE PRECISION A21, Q1
      DOUBLE PRECISION K1, W1, T1, X1

      DOUBLE PRECISION TEMP

      A21 = 1.0D+00
      Q1 = 1.0D+00

      H = (TN - T0) / DBLE(N)

      X(0) = X0
      T = T0

      TEMP = X(0)

      CALL TABLE_HEADER('RK-1 Time-Variant')
      CALL TABLE_ROW(0, T, TEMP)

      DO I = 1, N

            T = T0 + H * DBLE(I)

            T1 = T
            X1 = X(I-1)
            W1 = R8_NORMAL(SEED) * SQRT(Q1 * Q / H)
            K1 = H * F(X1, T1) + H * G(X1, T1) * W1
            X(I) = X1 + A21 * K1

            TEMP = X(I)

            CALL TABLE_ROW(I, T, TEMP)

      END DO

      CALL TABLE_FOOTER()

      RETURN
      END
