#include<param.h>
#include<poly_q.h>
#include<poly_red.h>
#include<poly_mult.h>
#include<sample.h>
#include<gaussian_avx.h>

void response(POLY_R *f,  POLY_R *f_J, POLY_R *z_a, POLY_R *z, POLY_R *z_J, POLY_R *a,  POLY_R  *b,  POLY_R *c,  POLY_R *d, uint64_t j_re[L], POLY_R *m, POLY_R *r, POLY_R *r_e,  POLY_R *r_f, POLY_R *r_g,  POLY_R *x){
    POLY_R *f_temp;
    POLY_R *z_a_temp;
    POLY_R *z_temp;
    POLY_R *f_J_temp;
    POLY_R *z_J_temp;

    // sample_alpha_prime(x);

    /* z_1 = y_1 + x * \rho */
    /* f = c + x * a*/
    for (int i = 0; i < 2 * L * NHAT * K; i++)
    {
        mult_r(f_temp + i, x, a + i);
        
        for (int j = 0; j < D; j++)
        {
            f[i].poly[j] = con_sub(f_temp[i].poly[j] + c[i].poly[j], Q);
        }
    }

    /* za = d + x * b*/
    for (int i = 0; i < M; i++)
    {
        mult_r(z_a_temp + i, x, b + i);
        
        for (int j = 0; j < D; j++)
        {
            z_a[i].poly[j] = con_sub(z_a_temp[i].poly[j] + d[i].poly[j], Q);
        }
    }

    /* z = rg + x * r*/
    for (int i = 0; i < M; i++)
    {
        mult_r(z_temp + i, x, r + i);
        
        for (int j = 0; j < D; j++)
        {
            z[i].poly[j] = con_sub(z_temp[i].poly[j] + r_g[i].poly[j], Q);
        }
    }

    //special
    /* fJ = m + x * j*/
    for (int i = 0; i < L; i++)
    {
        // mult_r(f_J_temp + i, x, j_re + i);
        for (int j = 0; j < D; j++)
        {
            f_J_temp[i].poly[j] = x->poly[j] * j_re[i];
            f_J[i].poly[j] = con_sub(f_J_temp[i].poly[j] + m[i].poly[j], Q);
        }
    }
    
    /* zJ = r_f + x * r_e*/
    for (int i = 0; i < M; i++)
    {
        mult_r(z_J_temp + i, x, r_e + i);
        
        for (int j = 0; j < D; j++)
        {
            z_J[i].poly[j] = con_sub(z_J_temp[i].poly[j] + r_f[i].poly[j], Q);
        }
    }
    

    //rej
    if(rej_op(z, z_a, z_J, z_temp, z_a_temp, z_J_temp) == 0){
        return 1;
    }

    if(rej_fJ(f_J, f_J_temp) == 0){
        return 1;
    }

    if(rej_f(f, f_temp) == 0){
        return 1;
    }

    

}   