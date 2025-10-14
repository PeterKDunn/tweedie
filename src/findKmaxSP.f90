
REAL(KIND=C_DOUBLE) FUNCTION findKmaxSP(i) BIND(C, NAME='findKmaxSP')
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN)    :: i
  REAL(KIND=C_DOUBLE)           :: tsmall, tlarge, abs1mp
  REAL(KIND=C_DOUBLE)           :: omegaInf, slope

  REAL(KIND=8)  :: current_y, current_mu, current_phi, pi

  ! --- Interface Declarations ---
  INTERFACE
    SUBROUTINE findImkd(i, t, Imdk)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
      INTEGER, INTENT(IN)               :: i
      REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imdk
    END SUBROUTINE findImkd
  END INTERFACE
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i


  ! Initialize
  abs1mp = ABS(1.0d00 - Cp)
  pi = 4.0D0 * DATAN(1.0D0)

  ! --- Core Logic ---
  
  ! We find a small-t approx, a large-t approx, and a combined approx.
  ! The SP should be the FIRST of these whose slope is NEGATIVE

  IF (CpSmall .AND. (current_y < current_mu) ) THEN
    ! Find the turning points.
    ! Try small-t approximation
    
    tsmall = SQRT(2.0d00 * (current_mu - current_y)/current_mu) * &
             current_mu**(1.0d00 - Cp) / current_phi
    findKmaxSP = tsmall
    
  ELSE
    omegaInf = (pi / 2.0d00) * &
               (1.0d00 - Cp)/(2.0d00*Cp - 1.0d00)
    tsmall = current_mu**(1.0d00 - Cp) / ( (1.0d00 - Cp)) * &
             TAN(omegaInf)
    
    ! --- Updated Call to findImkd (Crucial for stability) ---
    CALL findImkd(i, tsmall, slope)
    
    IF (slope <= 0.0d0) THEN
      findKmaxSP = tsmall
      RETURN
    END IF
  END IF
  
  ! Large-t approximation
  tlarge = (current_mu / current_y)**(Cp - 1.0d00) * &
           current_mu**(1.0d00 - Cp) / (current_phi * abs1mp)

  ! --- Updated Call to findImkd (Crucial for stability) ---
  CALL findImkd(i, tlarge, slope)

  IF (slope <= 0.0d0) THEN
    findKmaxSP = tlarge
    RETURN
  END IF

  ! Smooth interpolation: sum of small + large contributions
  findKmaxSP = tsmall + tlarge
  
END FUNCTION findKmaxSP
