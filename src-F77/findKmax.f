! This is a simple, module-less Fortran routine, mirroring the structure of findRek.
! The name MUST be in all caps (FINDKMAX) to guarantee the R linker finds 'findkmax_'.
SUBROUTINE FINDKMAX(p, mu, phi, k_max) 
    IMPLICIT NONE
    
    ! Arguments
    REAL(8), INTENT(IN) :: p, mu, phi
    INTEGER, INTENT(OUT) :: k_max  
    
    ! Local Variables
    REAL(8) :: exponent
    
    IF (mu <= 0.0D0) THEN
        k_max = 1
        RETURN
    END IF

    IF (p > 1.0D0 .AND. p < 2.0D0) THEN
        exponent = 1.0D0 / (p - 1.0D0)
        ! Calculate the result
        k_max = INT( 1.0D0 + (mu / phi) ** exponent )
        ! Ensure k_max is at least 10
        k_max = MAX(10, k_max + 1)
    ELSE
        k_max = 1
    END IF

END SUBROUTINE FINDKMAX
