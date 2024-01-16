        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 05 19:38:43 2020
        MODULE NEWT__genmod
          INTERFACE 
            FUNCTION NEWT(FUNC,X,FX,EPS,NITER) RESULT(IER)
              INTERFACE 
                FUNCTION FUNC(X,DFDX) RESULT(FX)
                  REAL(KIND=8) :: X(:)
                  REAL(KIND=8) ,OPTIONAL :: DFDX(:,:)
                  REAL(KIND=8) :: FX(SIZE(X))
                END FUNCTION FUNC
              END INTERFACE 
              REAL(KIND=8) :: X(:)
              REAL(KIND=8) :: FX(:)
              REAL(KIND=8), INTENT(IN) :: EPS
              INTEGER(KIND=4), INTENT(IN) :: NITER
              INTEGER(KIND=4) :: IER
            END FUNCTION NEWT
          END INTERFACE 
        END MODULE NEWT__genmod
