
MODULE TweedieIntZones
  
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL

  IMPLICIT NONE

CONTAINS 
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE integrateFirstRegion(i, mfirst, left_Of_Max, tmax, & ! INPUTS
                                    area0, zeroR)                          ! OUTPUT
    ! Integrates the initial region of the Fourier integral

    IMPLICIT NONE
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area0, zeroR
    INTEGER(C_INT), INTENT(IN)          :: mfirst, i
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax

    REAL(KIND=C_DOUBLE)                 :: t_Start_Point, zeroL, zeroBoundL, zeroBoundR
    REAL(KIND=C_DOUBLE)                 :: pi, TMP, current_y
    LOGICAL(C_BOOL)                     :: error


    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)

    ! Initialisation
    t_Start_Point = 0.0_C_DOUBLE
    zeroL = 0.0_C_DOUBLE
    zeroBoundL = 0.0_C_DOUBLE
    zeroBoundR = 0.0_C_DOUBLE

!WRITE(*,*) "left_Of_Max:", leftOfmax
    ! Find starting point for the first zero
    IF (left_Of_Max) THEN
      t_Start_Point = pi / current_y  
      zeroBoundL = 0_C_DOUBLE
      zeroBoundR = tmax   ! WAS: t_Start_Point * 2.0_C_DOUBLE
    ELSE
      ! Searching to the right of tmax
      t_Start_Point = tmax + pi / current_y  
      zeroBoundL = tmax
      zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
    END IF
!WRITE(*,*) "BOUNDS", zeroBoundL, t_Start_Point, zeroBoundR
    IF ( (t_Start_Point .GT. zeroBoundR) .OR. (t_Start_Point .LT. zeroBoundL) ) Then
      t_Start_Point = (zeroBoundL + zeroBoundR) / 2.0_c_DOUBLE
    END IF
  
    ! Find the zero
    zeroL = 0.0_C_DOUBLE
!WRITE(*,*) "BOUNDS---about to call find exact zeros", zeroBoundL, t_Start_Point, zeroBoundR

    CALL findExactZeros(i, mfirst, zeroBoundL, zeroBoundR, t_Start_Point, zeroR, left_Of_Max)
    ! findExactZeros may change the value of  left_Of_Max
!WRITE(*,*) "FOUND---zeroR", zeroR
    CALL evaluateImk(i, zeroR, TMP, error)
    IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, zeroR, 1)

    ! Find the area
    CALL GaussQuadrature(i, zeroL, zeroR, area0)
!WRITE(*,*) "From", zeroL, "to", zeroR, "(m=", m, "), area=", area0
    IF (Cverbose) THEN
      CALL DBLEPR("  *** INITIAL area:", -1, area0, 1 )
      CALL DBLEPR("         between t:", -1, 0.0_C_DOUBLE, 1 )
      CALL DBLEPR("             and t:", -1, zeroR, 1 )
      CALL INTPR( "   using (right) m:", -1, mfirst, 1)
    END IF

    END SUBROUTINE integrateFirstRegion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE integratePreAccRegions(i, m, mfirst, left_Of_Max, zeroL, tmax, mmax, flip_To_Other_Side, & ! INPUTS
                                   area1, zeroR, count_PreAcc_regions, converged_Pre)               ! OUTPUTS
    ! Integration of the region *after* the initial, but *before* acceleration is invoked.
    ! Potentially, everything converges in this step (without ever needing acceleration).
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area1
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroL, zeroR
    INTEGER(C_INT), INTENT(IN)          :: mfirst, mmax, i
    INTEGER(C_INT), INTENT(INOUT)       :: m
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max, flip_To_Other_Side
      ! flip_To_Other_Side signals that m has just crossed mmax;
      ! update integration side exactly once
    INTEGER(C_INT), INTENT(OUT)         :: count_PreAcc_Regions
    LOGICAL(C_BOOL), INTENT(OUT)        :: converged_Pre

    INTEGER(C_INT)                      :: mOld
    REAL(KIND=C_DOUBLE)                 :: zeroBoundL, zeroBoundR, zeroStartPoint, zero, sumA
    REAL(KIND=C_DOUBLE)                 :: area1Old, tolerance, sumAOld, current_y
    LOGICAL(C_BOOL)                     :: stop_PreAccelerate


    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)

    converged_Pre = .FALSE.           ! .TRUE. if it seems the integration has converged
    tolerance = 1.0E-12_C_DOUBLE      ! tolerance for concluding the area is not vhanging
    sumAOld = 0.0_C_DOUBLE            ! Previous estimation of the area
    count_PreAcc_Regions = 0_C_INT    ! Count how many pre-acc regions are evaluated
    area1 = 0.0_C_DOUBLE              ! Initalise
    mOld = m
    sumA = 0.0_C_DOUBLE
    
    IF (mfirst .EQ. -1_C_INT) THEN

      IF (Cverbose) CALL DBLEPR("  - No PRE-ACC area for y:", -1, current_y, 1 )
    ELSE ! That is,  mfirst  is not -1, so begins at zero *and* heads upwsards 
      CALL advanceM(i, m, mmax, mOld, left_Of_Max, flip_To_Other_Side)
      area1Old = 10.0_C_DOUBLE
      stop_PreAccelerate = .FALSE.

      DO WHILE ( .NOT.(stop_PreAccelerate ) )
        IF (left_Of_Max) THEN
          zeroBoundL = zeroR 
          zeroBoundR = tmax
          IF (flip_To_Other_Side) THEN 
            left_Of_Max = .FALSE.
            flip_To_Other_Side = .FALSE.
          END IF
        ELSE
          zeroBoundL = tmax
          zeroBoundR = zeroBoundL * 20.0_C_DOUBLE
        END IF
        zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
        zeroL = zeroR
        CALL improveKZeroBounds(i, m, left_Of_Max, zeroStartPoint, &
                                zeroBoundL, zeroBoundR)
        zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
        CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, &
                            zeroStartPoint, zero, left_Of_Max)
        zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
        zeroR = zero
  
        IF (count_PreAcc_Regions .GT. 1_C_INT) THEN
          area1Old = area1
          sumAOld = sumA
        END IF

        CALL GaussQuadrature(i, zeroL, zeroR, sumA)
!WRITE(*,*) "From", zeroL, "to", zeroR, "(m=", m, "), area=", sumA

        area1 = area1 + sumA
        count_PreAcc_Regions = count_PreAcc_Regions + 1

        IF (Cverbose) THEN 
          CALL DBLEPR("   Pre-acc:  Left t:", -1, zeroL, 1)
          CALL DBLEPR("   Pre-acc: Right t:", -1, zeroR, 1)
          CALL DBLEPR("               sumA:", -1, sumA, 1)
          CALL DBLEPR("              area1:", -1, area1, 1)
          CALL INTPR( "            with m::", -1, m, 1)
        END IF          
!CALL INTPR( "       >>>>>     with m::", -1, m, 1)
        
        ! Determine if we can stop pre-accelerating yet
        CALL stopPreAcc(i, tmax, zeroR, stop_PreAccelerate)

        mOld = m
        CALL advanceM(i, m, mmax, mOld, left_Of_Max, flip_To_Other_Side)
      END DO 
  
      IF (Cverbose) THEN
        CALL DBLEPR("  *** PRE-ACC area:", -1, sumA, 1 )
        CALL DBLEPR("         between t:", -1, zeroL, 1 )
        CALL DBLEPR("             and t:", -1, zeroR, 1 )
        CALL INTPR( "           using m:", -1, m, 1 )
      END IF
    END IF

  END SUBROUTINE integratePreAccRegions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE integrateAccelerationRegions(i, m, left_Of_Max, zeroL,  tmax, accMax, & ! INPUTS
                                          aimrerr, epsilon, xvec, wvec, Mmatrix, Nmatrix, & 
                                          mmax, min_Acc_Regions, exitstatus, flip_To_Other_Side, &
                                          convergence_Acc, &
                                          zeroR, converged_Accelerating, &
                                          relerr)
    ! Integration of the region *after* the initial and pre-acceleration regions, the tail area.
    ! Potentially, everything converges in pre-acceleration (without needing acceleration).


    IMPLICIT NONE

    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroR, relerr
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax, zeroL, aimrerr, epsilon
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: xvec(:), wvec(:), Mmatrix(:, :), Nmatrix(:, :)
    INTEGER(C_INT), INTENT(INOUT)       :: m, exitstatus
    INTEGER(C_INT), INTENT(IN)          :: mmax, min_Acc_Regions, accMax, i
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max, flip_To_Other_Side, convergence_Acc, converged_Accelerating

    INTEGER(C_INT)                      :: mOld, its_Acceleration
    REAL(KIND=C_DOUBLE)                 :: zeroBoundL, zeroBoundR, current_y
    REAL(KIND=C_DOUBLE)                 :: West, Wold, Wold2, pi, psi, zeroStartPoint
    LOGICAL(C_BOOL)                     :: keep_Accelerating


    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)

    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    
    ! Initialisation of acceleration 
    West = 3.0_C_DOUBLE                 ! The current estimate of the tail area
    Wold = 2.0_C_DOUBLE                 ! Wold and Wold2 are the previous two estimates of the tail area,
    Wold2 = 1.0_C_DOUBLE                !   to allow testing for convergence
    its_Acceleration = 0_C_INT          ! Number of integration regions used
    keep_Accelerating = .TRUE.          ! .TRUE. if we need more integration regions
    converged_Accelerating = .FALSE.    ! .TRUE. if it seems the integration has converged
    its_Acceleration = 0_C_INT          ! Number of acceleration integration regions
    
    ! This will be the very first, left-most value of t used in acceleration
    xvec(1) = zeroL

!      CALL advanceM(i, m, mmax, mOld, left_Of_Max, flip_To_Other_Side)
    mOld = m

    DO WHILE (keep_Accelerating )
      its_Acceleration = its_Acceleration + 1_C_INT
      mOld = m
      CALL advanceM(i, m, mmax, mOld, left_Of_Max, flip_To_Other_Side)
      zeroBoundL = zeroL
  
      IF (left_Of_Max) THEN
        zeroBoundR     = tmax
        zeroStartPoint = (zeroBoundL + zeroBoundR) / 2.0_C_DOUBLE
      ELSE
        ! Searching to the right of the maximum of Im k(t)
        IF (flip_To_Other_Side) THEN 
          ! Then this is the first time on the right side of the maximum
          left_Of_Max = .FALSE.
          flip_To_Other_Side = .FALSE.
          zeroBoundL = tmax
          zeroBoundR = tmax * 2.0_C_DOUBLE
          zeroStartPoint = (zeroBoundL + zeroBoundR) / 2.0_C_DOUBLE
        ELSE
          zeroStartPoint = zeroL + (pi / current_y)
          zeroBoundR = zeroStartPoint * 5.0_C_DOUBLE
        END IF
      END IF

      ! Find the exact zero
      CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zeroR, left_Of_Max)
      
      IF (its_Acceleration .GE. accMax) THEN
        CALL INTPR("Max acceleration regions reached. Stopping acceleration at m:", -1, m, 1)
        
        exitstatus = 1
        convergence_Acc = .TRUE.
        converged_Accelerating = .FALSE.
        EXIT
      END IF
  
      xvec(its_Acceleration + 1_C_INT) = zeroR
      
      CALL GaussQuadrature(i, zeroL, zeroR, psi)           ! psi: area of the latest region
!WRITE(*,*) "From", zeroL, "to", zeroR, "(m=", m, "), area=", psi
      wvec(its_Acceleration) = psi 
        ! wvec contains the sequence of integration areas, starting with the first and up to the limit
  
      ! Update past estimates
      Wold2 = Wold
      Wold  = West    
      ! Estimate the tail area using acceleration
      CALL accelerate(xvec, wvec, its_Acceleration, Mmatrix, Nmatrix, West)

      ! Check for convergence
      relerr = ( DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)
      IF (Cverbose) THEN
        CALL DBLEPR("  *** This acceleration area:", -1, psi, 1)
        CALL DBLEPR("                   between t:", -1, zeroL, 1)
        CALL DBLEPR("                       and t:", -1, zeroR, 1)
        CALL INTPR( "                     using m:", -1, m, 1)
        CALL DBLEPR("                with rel err:", -1, relerr, 1)
      END IF
  
      ! Declare convergence_Acc of we have sufficient regions, and relerr estimate is small
      IF ( (its_Acceleration .GE. min_Acc_Regions) .AND. &
           (relerr .LT. aimrerr) ) THEN
        keep_Accelerating = .FALSE.
        converged_Accelerating = .TRUE.
        convergence_Acc = .TRUE.
      END IF
      
      IF (its_Acceleration .GT. accMax) THEN
           keep_Accelerating = .FALSE.
           converged_Accelerating = .FALSE.
      END IF
!IF (m .LT. -13) STOP

!      zeroL = zeroR ! Next left zero, is previous region's right zero

    END DO  
  
  END SUBROUTINE integrateAccelerationRegions
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE stopPreAcc(i, tmax, zeroL, stop_PreAccelerate)
    ! Determine if it is OK to stop pre-accelerating, and start using acceleration
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN) :: tmax
    REAL(KIND=C_DOUBLE), INTENT(IN) :: zeroL
    INTEGER (C_INT), INTENT(IN)     :: i
    LOGICAL(C_BOOL), INTENT(OUT)    :: stop_PreAccelerate
    
    INTEGER (C_INT)       :: nmax
    REAL(KIND=C_DOUBLE)   :: MM, Rek, Rekd, tstop, pi
    REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi
    
    
    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    
    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)
    current_mu   = Cmu(i)
    current_phi  = Cphi(i)


    ! Stop condition for pre-acceleration.
    ! Ensure that we have passed the peak of Im k(t), so that acceleration can be used
    stop_PreAccelerate = .FALSE.
    IF ( Cp .GT. 2.0_C_DOUBLE) THEN
      IF (current_y .GT. current_mu) THEN
        stop_PreAccelerate = .TRUE.
      ELSE
        IF (zeroL .GT. tmax) THEN
          stop_PreAccelerate = .TRUE.
        END IF
      END IF
    ELSE
      IF (current_y .GT. current_mu) THEN
        stop_PreAccelerate = .TRUE.
      ELSE
        MM = 1.0_C_DOUBLE / (2.0_C_DOUBLE  * (Cp - 1.0_C_DOUBLE))
        IF ( MM .LE. 1.0_C_DOUBLE ) THEN
          stop_PreAccelerate = .TRUE.
        ELSE
        ! Check when t is larger than the last turning point of Re k(t) 
        ! Or when exp{Re(k)/t} is so small that it makes no difference... but 
        ! care is needed: Re k(t) is not necessarily convex here
           IF ( ABS( MM - FLOOR(MM) ) < 1.0E-09_C_DOUBLE ) then
              nmax = FLOOR(MM) - 1_C_INT
           ELSE
              nmax = FLOOR(MM)
           END IF
           tstop = current_mu**(1.0_C_DOUBLE - Cp) / ((1.0_C_DOUBLE - Cp) * current_phi) *   & 
                   DTAN( DBLE(nmax) * pi * (1.0_C_DOUBLE - Cp) )
           IF (zeroL .GT. tstop) stop_PreAccelerate = .TRUE.
        END IF
      END IF
    END IF
    
    ! Sometimes this takes forever to flag stop_preacc as TRUE,
    ! so also check if exp{Re k(t)/t} is small
    CALL evaluateRek( i, zeroL, Rek)
    CALL evaluateRekd(i, zeroL, Rekd)
    
    IF (zeroL .GT. 0.0_C_DOUBLE) THEN
      IF ( ( (DEXP(Rek)/zeroL) .LT. 1.0E-07_C_DOUBLE) .AND.          & 
          (Rekd .LT. 0.0_C_DOUBLE) ) then
       stop_PreAccelerate = .TRUE.
        IF (Cverbose) THEN 
          CALL DBLEPR("Pre-accelerating stopping. Rek(t) small:", -1, (DEXP(Rek)/zeroL), 1)
        END IF
      END IF
    END IF

    IF (zeroL .GT. 0.0_C_DOUBLE) THEN
      IF ( ( (DEXP(Rek)/zeroL) .LT. 1.0E-15_C_DOUBLE)  ) then
        stop_PreAccelerate = .TRUE.
        IF (Cverbose) THEN 
          CALL DBLEPR("Pre-accelerating stopping. Rek(t) small:", -1, (DEXP(Rek)/zeroL), 1)
        END IF
      END IF
    END IF

  END SUBROUTINE stopPreAcc
    
END MODULE TweedieIntZones
    
