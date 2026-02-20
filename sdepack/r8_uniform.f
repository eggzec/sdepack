      FUNCTION R8_UNIFORM(seed)
C     *******************************************************************
C
C       R8_UNIFORM generates a pseudorandom double-precision value
C       uniformly distributed in the interval (0, 1).
C
C       Discussion:
C
C         This routine implements the Park-Miller linear congruential
C         generator (Lehmer generator):
C
C           seed = 16807 * seed mod (2**31 - 1)
C           r8_uniform = seed / (2**31 - 1)
C
C         The integer arithmetic is done using Schrage's method to
C         avoid overflow.
C
C       Licensing:
C
C         This code is distributed under the MIT license.
C
C       Modified:
C
C         17 July 2006
C
C       Author:
C
C         John Burkardt
C
C       References:
C
C         Paul Bratley, Bennett Fox, Linus Schrage,
C         A Guide to Simulation,
C         Springer Verlag, pages 201-202, 1983.
C
C         Pierre L'Ecuyer,
C         Random Number Generation,
C         in Handbook of Simulation,
C         edited by Jerry Banks,
C         Wiley Interscience, page 95, 1998.
C
C       Parameters:
C
C         Input/output:
C           SEED   - INTEGER
C                    On input, a seed value. Must NOT be 0.
C                    On output, the updated seed value.
C
C         Output:
C           R8_UNIFORM_01 - DOUBLE PRECISION
C                    A pseudorandom value uniformly distributed
C                    in the interval (0, 1).
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER k
      DOUBLE PRECISION r8_uniform
      INTEGER seed

      k = seed / 127773

      seed = 16807 * (seed - k * 127773) - k * 2836

      IF (seed .LT. 0) THEN
        seed = seed + 2147483647
      END IF

C
C     Although SEED can be represented exactly as a 32-bit integer,
C     it generally cannot be represented exactly as a 32-bit real.
C
      r8_uniform = DBLE(seed) * 4.656612875D-10

      RETURN
      END
