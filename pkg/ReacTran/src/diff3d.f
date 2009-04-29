C bind two matrices to an array in x-direction; mat1: first, mat2:last

        SUBROUTINE mbindx2(Nx,Ny,Nz,ARRAY,MAT1,MAT2,ARR)

        INTEGER Nx, Ny, Nz
        DOUBLE PRECISION ARRAY(Nx,Ny,Nz), ARR (Nx+2,Ny,Nz)
        DOUBLE PRECISION MAT1(Ny,Nz), MAT2(Ny,Nz)

        INTEGER I, J, K
        
         DO J = 1,Ny
           DO K = 1, Nz
             ARR(1,J,K)    = MAT1(J,K)
             ARR(Nx+2,J,K) = MAT2(J,K)
             DO I = 1, Nx
              ARR(I+1,J,K) = ARRAY(I,J,K)
             ENDDO
           ENDDO
         ENDDO

        RETURN
        END

C bind two matrices to an array in Y-direction; mat1: first, mat2:last

        SUBROUTINE mbindy2(Nx,Ny,Nz,ARRAY,MAT1,MAT2,ARR)

        INTEGER Nx, Ny, Nz
        DOUBLE PRECISION ARRAY(Nx,Ny,Nz), ARR (Nx,Ny+2,Nz)
        DOUBLE PRECISION MAT1(Nx,Nz), MAT2(Nx,Nz)

        INTEGER I, J, K

         DO I = 1,Nx
           DO K = 1, Nz
             ARR(I,1,K)    = MAT1(I,K)
             ARR(I,Ny+2,K) = MAT2(I,K)
             DO J = 1, Ny
              ARR(I,J+1,K) = ARRAY(I,J,K)
             ENDDO
           ENDDO
         ENDDO

        RETURN
        END

C bind two matrices to an array in Z-direction; mat1: first, mat2:last

        SUBROUTINE mbindz2(Nx,Ny,Nz,ARRAY,MAT1,MAT2,ARR)

        INTEGER Nx, Ny, Nz
        DOUBLE PRECISION ARRAY(Nx,Ny,Nz), ARR (Nx,Ny,Nz+2)
        DOUBLE PRECISION MAT1(Nx,Ny), MAT2(Nx,Ny)

        INTEGER I, J, K

         DO I = 1,Nx
           DO J = 1, Ny
             ARR(I,J,1)    = MAT1(I,J)
             ARR(I,J,Nz+2) = MAT2(I,J)
             DO K = 1, Nz
              ARR(I,J,K+1) = ARRAY(I,J,K)
             ENDDO
           ENDDO
         ENDDO

        RETURN
        END



C bind one matrix to an array in x-direction; mat

        SUBROUTINE mbindx1(Nx,Ny,Nz,ARRAY,MAT,ARR,Typ)

        INTEGER Nx, Ny, Nz, Typ
        DOUBLE PRECISION ARRAY(Nx,Ny,Nz), ARR (Nx+1,Ny,Nz)
        DOUBLE PRECISION MAT1(Ny,Nz), MAT2(Ny,Nz)

        INTEGER I, J, K

         DO J = 1,Ny
           DO K = 1, Nz
             IF (Typ==1) THEN
               ARR(1,J,K)    = MAT(J,K)
             ELSE
               ARR(Nx+1,J,K) = MAT(J,K)
             ENDIF
             DO I = 1, Nx
              ARR(I+1,J,K) = ARRAY(I,J,K)
             ENDDO
           ENDDO
         ENDDO

        RETURN
        END

C bind one matrix to an array in Y-direction; mat

        SUBROUTINE mbindy1(Nx,Ny,Nz,ARRAY,MAT,ARR,Typ)

        INTEGER Nx, Ny, Nz, Typ
        DOUBLE PRECISION ARRAY(Nx,Ny,Nz), ARR (Nx,Ny+1,Nz)
        DOUBLE PRECISION MAT1(Nx,Nz), MAT2(Nx,Nz)

        INTEGER I, J, K

         DO I = 1,Nx
           DO K = 1, Nz
             IF (typ == 1) THEN
               ARR(I,1,K)    = MAT(I,K)
             ELSE
               ARR(I,Ny+1,K) = MAT(I,K)
             ENDIF
             DO J = 1, Ny
              ARR(I,J+1,K) = ARRAY(I,J,K)
             ENDDO
           ENDDO
         ENDDO

        RETURN
        END

C bind two matrices to an array in Z-direction; mat1: first, mat2:last

        SUBROUTINE mbindz1(Nx,Ny,Nz,ARRAY,MAT,ARR,Typ)

        INTEGER Nx, Ny, Nz, Typ
        DOUBLE PRECISION ARRAY(Nx,Ny,Nz), ARR (Nx,Ny,Nz+1)
        DOUBLE PRECISION MAT1(Nx,Ny), MAT2(Nx,Ny)

        INTEGER I, J, K

         DO I = 1,Nx
           DO J = 1, Ny
             IF (typ == 1) THEN
               ARR(I,J,1)    = MAT(I,J)
             ELSE
               ARR(I,J,Nz+1) = MAT(I,J)
             ENDIF
             DO K = 1, Nz
              ARR(I,J,K+1) = ARRAY(I,J,K)
             ENDDO
           ENDDO
         ENDDO

        RETURN
        END
