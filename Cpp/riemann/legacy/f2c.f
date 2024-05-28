! f2c.f
! Example in fixed source format.
      PROGRAM F2C
      IMPLICIT NONE
      REAL :: C, F

      PRINT *, '  Fahrenheit     Celsius'
      PRINT *, '--------------------------'

! Output table:
      DO F = 30, 220, 10
      C = (5.0 / 9.0) * (F - 32.0)
      PRINT '(F13.1, F12.3)', F, C
      END DO
      END PROGRAM F2C
