      SUBROUTINE TABLE_HEADER(TITLE)
C     *****************************************************************
C
C       TABLE_HEADER prints a formatted table header for SDE output.
C
C     *****************************************************************

      IMPLICIT NONE

      CHARACTER*(*) TITLE

      WRITE(*,'(A)') ' '
      WRITE(*,'(A)') ' ' // TITLE
      WRITE(*,'(A)')
     &  ' +----------+----------------+----------------+'
      WRITE(*,'(A)')
     &  ' |   STEP   |        T       |        X       |'
      WRITE(*,'(A)')
     &  ' +----------+----------------+----------------+'

      RETURN
      END
      SUBROUTINE TABLE_ROW(I, T, X)
C     *****************************************************************
C
C       TABLE_ROW prints one row in the SDE solution table.
C
C     *****************************************************************

      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION T
      DOUBLE PRECISION X

      WRITE(*,'(A,I8,A,F14.6,A,G14.6,A)')
     &  ' | ', I, ' | ', T, ' | ', X, ' |'

      RETURN
      END
      SUBROUTINE TABLE_FOOTER()
C     *****************************************************************
C
C       TABLE_FOOTER prints a table footer for SDE output.
C
C     *****************************************************************

      IMPLICIT NONE

      WRITE(*,'(A)')
     &  ' +----------+----------------+----------------+'

      RETURN
      END
