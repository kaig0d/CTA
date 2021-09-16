#ifndef MERKLE_TREE_VERIFY_H


#include<param.h>
#include<poly_q.h>

int verify(POLY_QHAT *p, POLY_Q *C_upp, POLY_Q *D_upp, POLY_Q *E_upp, POLY_Q *F_upp, POLY_Q *G_0, POLY_Q *G_1, POLY_Q *ConstTPrime, POLY_Q *W_upp, POLY_Q *x, POLY_Q *y, POLY_Q * gamma, POLY_Q *f,  POLY_Q *f_J, POLY_Q *z_a, POLY_Q *z, POLY_Q *z_J, POLY_Q mat1[][N], POLY_Q mat2[][1], POLY_Q H_0_HAT[][L * NHAT], POLY_Q H_1_HAT[][L * NHAT], POLY_Q NHAT_mat[][L * NHAT]);

#endif