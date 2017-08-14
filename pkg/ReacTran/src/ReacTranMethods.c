#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif
/*     subroutine advection(N, Y, dt, h, hint, v, Bcup, Bcdown,              &
     &       Yup, Ydown, VFint, VF, Aint, A, method, mode, split,            &
     &       dY, cu, it 
     
       subroutine advectvol(N, Y, dt, V, Vint, flow, Bcup, Bcdown,           &
     &                     Yup, Ydown,  method,mode,dY, cu, it)
     */

void F77_NAME(advection)(int*, double*, double*, double*, double*, double *, 
                         int*, int*, 
                         double*, double*, double*, double*, double*, double*, 
                         int*, int*, int*, double*, double*, int*);

void F77_NAME(advectvol)(int*, double*, double*, double*, double*, double *, 
                         int*, int*, 
                         double*, double*, int*, int*, double*, double*, int*);
        
R_FortranMethodDef ReacfortranMethods[] = {
 {"advection",   (DL_FUNC) &F77_SUB(advection), 20},
 {"advectvol",   (DL_FUNC) &F77_SUB(advectvol), 15},
 {NULL, NULL, 0}
};

void R_init_ReacTran(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, ReacfortranMethods, NULL);
}
