        !COMPILER-GENERATED INTERFACE MODULE: Tue Aug 11 16:47:04 2020
        MODULE CUADRATICA__genmod
          INTERFACE 
            FUNCTION CUADRATICA(X,DFDX) RESULT(FUNX)
              REAL(KIND=8) :: X(:)
              REAL(KIND=8) ,OPTIONAL :: DFDX(:,:)
              REAL(KIND=8) :: FUNX(SIZE(X))
            END FUNCTION CUADRATICA
          END INTERFACE 
        END MODULE CUADRATICA__genmod