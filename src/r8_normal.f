      FUNCTION R8_NORMAL(SEED)
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
C         the code saVes the second value for the next call.
C
C         If the user changes the SEED between calls, the routine
C         automatically resets itself and discards the saVed data.
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
C                    SEED for the random number generator.
C                    Updated internally.
C
C         Output:
C           R8_NORMAL - DOUBLE PRECISION
C                    A pseudorandom variate from the standard normal
C                    distribution N(0,1).
C
C     *******************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.141592653589793D+00)
      DOUBLE PRECISION R1
      DOUBLE PRECISION R2
      DOUBLE PRECISION R8_NORMAL
      DOUBLE PRECISION R8_UNIFORM
      INTEGER SEED
      INTEGER SEED1
      INTEGER SEED2
      INTEGER SEED3
      INTEGER USED
      DOUBLE PRECISION V1
      DOUBLE PRECISION V2

      SAVE SEED1
      SAVE SEED2
      SAVE SEED3
      SAVE USED
      SAVE V2

      DATA SEED2 / 0 /
      DATA USED / 0 /
      DATA V2 / 0.0D+00 /

C
C     If USED is odd, but the input SEED does not match
C     the output SEED on the previous call, then the user has changed
C     the SEED.  Wipe out internal memory.
C
      IF (MOD(USED, 2) .EQ. 1) THEN

        IF (SEED .NE. SEED2) THEN
          USED = 0
          SEED1 = 0
          SEED2 = 0
          SEED3 = 0
          V2 = 0.0D+00
        END IF

      END IF

C
C     If USED is even, generate two uniforms, create two normals,
C     return the first normal and its corresponding SEED.
C
      IF (MOD(USED, 2) .EQ. 0) THEN

        SEED1 = SEED

        R1 = R8_UNIFORM(SEED)

        IF (R1 .EQ. 0.0D+00) THEN
          WRITE (*,'(a)') ' '
          WRITE (*,'(a)') 'R8_NORMAL - Fatal error!'
          WRITE (*,'(a)') '  R8_UNIFORM returned a value of 0.'
          STOP
        END IF

        SEED2 = SEED

        R2 = R8_UNIFORM(SEED)

        SEED3 = SEED

        V1 = SQRT(-2.0D+00 * LOG(R1)) * COS(2.0D+00 * PI * R2)
        V2 = SQRT(-2.0D+00 * LOG(R1)) * SIN(2.0D+00 * PI * R2)

        R8_NORMAL = V1
        SEED = SEED2

C
C     If USED is odd (and the input SEED matched the output value from
C     the previous call), return the second normal and its corresponding
C     SEED.
C
      ELSE

        R8_NORMAL = V2
        SEED = SEED3

      END IF

      USED = USED + 1

      RETURN
      END
