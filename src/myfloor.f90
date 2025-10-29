FUNCTION myfloor(x) RESULT(r_myfloor) 
  ! A floor function, for my purposes.
  USE ISO_C_BINDING, ONLY: C_DOUBLE
  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: x
  INTEGER                           :: r_myfloor


  IF (x .GE. 0.0_C_DOUBLE) THEN
    ! For positive numbers (2.5 -> 2, 2.0 -> 2)
    r_myfloor = INT(x)
  ELSE
    ! For negative numbers:
    ! We must check if it's already an integer before subtracting 1
    IF (x .EQ. REAL(INT(x), KIND=C_DOUBLE)) THEN
        r_myfloor = INT(x)  ! If x is -2.0, result is -2
    ELSE
        r_myfloor = INT(x) - 1 ! If x is -2.5, result is -3
    END IF
  END IF

END FUNCTION myfloor
