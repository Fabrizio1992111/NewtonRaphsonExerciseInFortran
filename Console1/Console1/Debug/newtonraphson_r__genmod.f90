        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 05 19:13:53 2020
        MODULE NEWTONRAPHSON_R__genmod
          INTERFACE 
            FUNCTION NEWTONRAPHSON_R(FUNC,X,FX,EPS,NITER) RESULT(IER)
              INTERFACE 
                FUNCTION FUNC(X,DFDX) RESULT(FX)
                  REAL(KIND=4) :: X(:)
                  REAL(KIND=4) ,OPTIONAL :: DFDX(:,:)
                  REAL(KIND=4) :: FX(SIZE(X))
                END FUNCTION FUNC
              END INTERFACE 
              REAL(KIND=4), INTENT(INOUT) :: X(:)
              REAL(KIND=4) ,OPTIONAL, INTENT(OUT) :: FX(:)
              REAL(KIND=4), INTENT(IN) :: EPS
              INTEGER(KIND=4), INTENT(IN) :: NITER
              INTEGER(KIND=4) :: IER
            END FUNCTION NEWTONRAPHSON_R
          END INTERFACE 
        END MODULE NEWTONRAPHSON_R__genmod
