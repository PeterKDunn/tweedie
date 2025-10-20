
SUBROUTINE findAccelStart(i, tRightMost) 
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: tRightMost      ! The output starting point for acceleration
  
  REAL(KIND=C_DOUBLE)              :: pi, omegaRM, lRightMostD
  INTEGER                   :: lRightMost
  INTEGER(C_INT)            :: i
  REAL(KIND=C_DOUBLE)       :: current_y, current_mu, current_phi
  
  INTEGER, EXTERNAL :: myfloor

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

  IF ( (Cp > 1.4_C_DOUBLE) .AND. (Cp < 1.6_C_DOUBLE) ) THEN
    tRightMost = 1.0_C_DOUBLE
  ELSE
    lRightMostD = DACOS( (current_y / current_mu)**(Cp - 1.0_C_DOUBLE) ) / (pi * (Cp - 1.0_C_DOUBLE) )
                  
    lRightMost = myfloor(lRightMostD) + 1

    ! Final calculation for tRightMost
    omegaRM = REAL(lRightMost, KIND=C_DOUBLE) * pi * (1.0_C_DOUBLE - Cp)
    tRightMost = DTAN(omegaRM) * current_mu ** (1.0_C_DOUBLE - Cp) / &
                  ( current_phi * (1.0_C_DOUBLE - Cp) )
  END IF
  
END SUBROUTINE findAccelStart
