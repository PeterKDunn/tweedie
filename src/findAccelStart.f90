
SUBROUTINE findAccelStart(i, tRightMost) 
  USE tweedie_params_mod

  IMPLICIT NONE
  
  REAL(KIND=8), INTENT(OUT) :: tRightMost      ! The output starting point for acceleration
  
  REAL(KIND=8)              :: pi, omegaRM, lRightMostD
  INTEGER                   :: lRightMost, i
  REAL(KIND=8)              :: current_y, current_mu, current_phi
  
  INTEGER, EXTERNAL :: myfloor

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  pi = 4.0D0 * DATAN(1.0D0)

  IF ((Cp > 1.4d00) .AND. (Cp < 1.6d00)) THEN
    tRightMost = 1.0d0
  ELSE
    lRightMostD = ACOS( (current_y / current_mu)**(Cp - 1.0d0) ) / &
                  (pi * (Cp - 1.0d0) )
                  
    lRightMost = myfloor(lRightMostD) + 1

    ! Final calculation for tRightMost
    omegaRM = REAL(lRightMost, KIND=8) * pi * (1.0d0 - Cp)
    tRightMost = TAN(omegaRM) * current_mu ** (1.0d0 - Cp) / &
                 ( current_phi * (1.0d0 - Cp) )
  END IF
  
END SUBROUTINE findAccelStart
