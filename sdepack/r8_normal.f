      FUNCTION R8_NORMAL(seed)
C     *******************************************************************
C
C       R8_NORMAL generates a unit pseudonormal double-precision
C       value using the Box-Muller method.
C
C       Discussion:
C
C         This routine uses the Box Muller method, which requires pairs
C         of uniform random values to generate a pair of normal random
C         values. Since the user typically requests values one at a time,
C         the code saves the second value for the next call.
C
C         If the user changes the SEED between calls, the routine
C         automatically resets itself and discards the saved data.
C
C       Licensing:
C
C         This code is distributed under the MIT license.
C
C       Modified:
C
C         08 January 2007
C
C       Author:
C
C         John Burkardt
C
C       Parameters:
C
C         Input/output:
C           SEED   - INTEGER
C                    Seed for the random number generator.
C                    Updated internally.
C
C         Output:
C           R8_NORMAL_01 - DOUBLE PRECISION
C                    A pseudorandom variate from the standard normal
C                    distribution N(0,1).
C
C     *******************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION pi
      PARAMETER (pi = 3.141592653589793D+00)
      DOUBLE PRECISION r1
      DOUBLE PRECISION r2
      DOUBLE PRECISION r8_normal
      DOUBLE PRECISION r8_uniform
      INTEGER seed
      INTEGER seed1
      INTEGER seed2
      INTEGER seed3
      INTEGER used
      DOUBLE PRECISION v1
      DOUBLE PRECISION v2

      SAVE seed1
      SAVE seed2
      SAVE seed3
      SAVE used
      SAVE v2

      DATA seed2 / 0 /
      DATA used / 0 /
      DATA v2 / 0.0D+00 /

C
C     If USED is odd, but the input SEED does not match
C     the output SEED on the previous call, then the user has changed
C     the seed.  Wipe out internal memory.
C
      IF (MOD(used, 2) .EQ. 1) THEN

        IF (seed .NE. seed2) THEN
          used = 0
          seed1 = 0
          seed2 = 0
          seed3 = 0
          v2 = 0.0D+00
        END IF

      END IF

C
C     If USED is even, generate two uniforms, create two normals,
C     return the first normal and its corresponding seed.
C
      IF (MOD(used, 2) .EQ. 0) THEN

        seed1 = seed

        r1 = r8_uniform(seed)

        IF (r1 .EQ. 0.0D+00) THEN
          WRITE (*,'(a)') ' '
          WRITE (*,'(a)') 'R8_NORMAL - Fatal error!'
          WRITE (*,'(a)') '  R8_UNIFORM returned a value of 0.'
          STOP
        END IF

        seed2 = seed

        r2 = r8_uniform(seed)

        seed3 = seed

        v1 = SQRT(-2.0D+00 * LOG(r1)) * COS(2.0D+00 * pi * r2)
        v2 = SQRT(-2.0D+00 * LOG(r1)) * SIN(2.0D+00 * pi * r2)

        r8_normal = v1
        seed = seed2

C
C     If USED is odd (and the input SEED matched the output value from
C     the previous call), return the second normal and its corresponding
C     seed.
C
      ELSE

        r8_normal = v2
        seed = seed3

      END IF

      used = used + 1

      RETURN
      END
