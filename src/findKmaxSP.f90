
REAL(KIND=C_DOUBLE) FUNCTION findKmaxSP(i) 
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)    :: i
  REAL(KIND=C_DOUBLE)           :: tsmall, tlarge, abs1mp
  REAL(KIND=C_DOUBLE)           :: omegaInf, slope, pi

  REAL(KIND=C_DOUBLE)           :: current_y, current_mu, current_phi

  INTERFACE
    SUBROUTINE findImkd(i, t, Imkd)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
      INTEGER(C_INT), INTENT(IN)         :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)   :: Imkd
    END SUBROUTINE findImkd
  END INTERFACE
  

  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i

  ! Initialize
  abs1mp = ABS(1.0_C_DOUBLE - Cp)
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

  ! We find a small-t approx, a large-t approx, and a combined approx.
  ! The SP should be the FIRST of these whose slope is NEGATIVE

  IF (CpSmall .AND. (current_y < current_mu) ) THEN
    ! CASE: 1 < p < 2, AND y < mu

    ! Try small-t approximation
    tsmall = DSQRT(2.0_C_DOUBLE * (current_mu - current_y)/current_mu) * &
             current_mu**(1.0_C_DOUBLE - Cp) / current_phi
    findKmaxSP = tsmall
  ELSE
    ! This should be the CASE: p > 2, AND y < mu
    omegaInf = (pi / 2.0_C_DOUBLE) * &
               (1.0_C_DOUBLE - Cp)/(2.0_C_DOUBLE*Cp - 1.0_C_DOUBLE)
    tsmall = current_mu**(1.0_C_DOUBLE - Cp) / ( (1.0_C_DOUBLE - Cp)) * &
             DTAN(omegaInf)
    
    CALL findImkd(i, tsmall, slope)
    
    IF (slope .LE. 0.0_C_DOUBLE) THEN
      findKmaxSP = tsmall
      RETURN
    END IF
  END IF
  
  ! Large-t approximation
  tlarge = (current_mu / current_y)**(Cp - 1.0_C_DOUBLE) * &
            current_mu**(1.0_C_DOUBLE - Cp) / (current_phi * abs1mp)

  CALL findImkd(i, tlarge, slope)

  IF (slope .LE. 0.0_C_DOUBLE) THEN
    findKmaxSP = tlarge
    RETURN
  END IF

  ! Sum of small + large contributions
  findKmaxSP = tsmall + tlarge
  
END FUNCTION findKmaxSP
