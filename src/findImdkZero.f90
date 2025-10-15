
SUBROUTINE findImdkZero(i, t, f, df) 
  USE tweedie_params_mod, ONLY: Cphi, Cmu, Cy

  IMPLICIT NONE
  
  ! We must declare the interfaces for the subroutines we call, 
  ! as we are not using MODULEs. This is F90 best practice for clarity.
  INTERFACE
    SUBROUTINE findImkd(i, t, Imdk)

      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)   :: t
      INTEGER, INTENT(IN)        :: i
      REAL(KIND=8), INTENT(OUT)  :: Imdk
    END SUBROUTINE findImkd
    
    SUBROUTINE findImkdd(i, t, Imddk)

      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)   :: t
      INTEGER, INTENT(IN)        :: i
      REAL(KIND=8), INTENT(OUT)  :: Imddk
    END SUBROUTINE findImkdd
  END INTERFACE
  
  REAL(KIND=8), INTENT(IN)     :: t
  REAL(KIND=8), INTENT(OUT)    :: f, df
  INTEGER, INTENT(IN)          :: i

  
  ! Local Variables
  REAL(KIND=8)  :: Imdk, Imddk
  REAL(KIND=8)  :: current_y, current_mu, current_phi
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i


  ! Call to findImkd is now correct (6 arguments)
  CALL findImkd(i, t, Imdk)
  
  ! Call to findImkdd is now correct (6 arguments)
  CALL findImkdd(i, t, Imddk)

  ! Assign results
  f  = Imdk
  df = Imddk

END SUBROUTINE findImdkZero
