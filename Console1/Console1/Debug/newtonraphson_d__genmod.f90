        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 05 23:14:49 2020
        MODULE NEWTONRAPHSON_D__genmod
          INTERFACE 
            FUNCTION NEWTONRAPHSON_D(FUNC,X,FX,EPS,NITER) RESULT(IER)
              INTERFACE 
                FUNCTION FUNC(X,DFDX) RESULT(FX)
                  REAL(KIND=8) :: X(:)
                  REAL(KIND=8) ,OPTIONAL :: DFDX(:,:)
                  REAL(KIND=8) :: FX(SIZE(X))
                END FUNCTION FUNC
              END INTERFACE 
              REAL(KIND=8), INTENT(INOUT) :: X(:)
              REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: FX(:)
              REAL(KIND=8), INTENT(IN) :: EPS
              INTEGER(KIND=4), INTENT(IN) :: NITER
              INTEGER(KIND=4) :: IER
            END FUNCTION NEWTONRAPHSON_D
          END INTERFACE 
        END MODULE NEWTONRAPHSON_D__genmod
