#include<stdio.h>
#include<stdint.h>
#include<sammat.h>
#include<sample.h>
#include<param.h>
#include<poly_q.h>
#include<merkle_tree_commit.h>
#include<merkle_tree_response.h>
#include<merkle_tree_verify.h>


int main(){
    POLY_Q *b;
    POLY_Q *r_e;
    POLY_Q *r;
    POLY_Q *c;
    POLY_Q *m;
    POLY_Q *d;
    POLY_Q *r_f;
    POLY_Q *r_g;
    POLY_Q *C_upp;
    POLY_Q *D_upp;
    POLY_Q *E_upp;
    POLY_Q *F_upp;
    POLY_Q *G_0;
    POLY_Q *W_upp;
    POLY_Q *G_1;
    POLY_Q *ConstTPrime;
    //
    POLY_Q *a;
    uint64_t j_re[L];
    POLY_Q *p;
    POLY_Q  *y;
    POLY_R *x;
    POLY_Q *gamma;
    POLY_Q mat1[2 * L * NHAT * K][N];
    POLY_Q mat2[M][1];
    POLY_Q H_0[NHAT * K][NHAT];
    POLY_Q H_1[NHAT * K][NHAT];
    POLY_Q H_0_HAT[2 * L * NHAT * K][L * NHAT];
    POLY_Q H_1_HAT[2 * L * NHAT * K][L * NHAT];
    POLY_Q N_mat[NHAT * K][NHAT];
    POLY_Q IN_mat[2 * (L - 1) * NHAT *K][(L - 1) * NHAT];
    POLY_Q NHAT_mat[2 * L * NHAT * K][L * NHAT];
    
    POLY_R *f;
    POLY_R *f_J;
    POLY_R *z_a;
    POLY_R *z;
    POLY_R *z_J;
    
    sample_alpha_L(gamma);
    sample_alpha_NL(y);
    sample_alpha_prime(x);
    sample_mat1(mat1);
    sample_mat2(mat2);
    sample_mat3(H_0);
    sample_mat3(H_1);

    for(int i = 0 ; i < L; i++){
        for(int j = i * 2 * NHAT * K; j < (2 * i + 1)* NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                for(int l = 0; l < D; l++){
                    H_0_HAT[j][k].poly[l] = H_0[j - i * 2 * NHAT * K][k - i * NHAT].poly[l];

                }

            }
        }
    }

    for(int i = 0 ; i < L; i++){
        for(int j = (2 * i + 1)* NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                for(int l = 0; l < D; l++){
                    H_0_HAT[j][k].poly[l] = H_1[j - (2 * i + 1)* NHAT * K][k - i * NHAT].poly[l];
                    
                }

            }
        }
    }


    for(int i = 0 ; i < L; i++){
        for(int j = i * 2 * NHAT * K; j < (2 * i + 1)* NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                for(int l = 0; l < D; l++){
                    H_1_HAT[j][k].poly[l] = H_1[j - i * 2 * NHAT * K][k - i * NHAT].poly[l];

                }

            }
        }
    }

    for(int i = 0 ; i < L; i++){
        for(int j = (2 * i + 1)* NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                for(int l = 0; l < D; l++){
                    H_1_HAT[j][k].poly[l] = H_0[j - (2 * i + 1)* NHAT * K][k - i * NHAT].poly[l];
                    
                }

            }
        }
    }

    for(int i = 0; i < NHAT; i++){
        for(int j = i * K; j < (i + 1) * K; j = j++){
            N_mat[j][i].poly[0] = 1 << (j - i * K);

        }
        
    }


    for(int i = 0 ; i < L - 1; i++){
        for(int j = i * 2 * NHAT * K; j < (2 * i + 1)* NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                for(int l = 0; l < D; l++){
                    IN_mat[j][k].poly[l] = N_mat[j - i * 2 * NHAT * K][k - i * NHAT].poly[l];

                }

            }
        }
    }

    for(int i = 2 * NHAT * K; i < 2 * L * NHAT * K; i++){
        for(int j = 0; j < (L - 1) * NHAT; j++){
            for(int k = 0; k < D; k++){
                NHAT_mat[i][j].poly[k] = IN_mat[i - 2 * NHAT * K][j].poly[k];

            }
        }
    }


    //just sample, not right!!!
    sample_1(a, 2 * L * NHAT * K);
    sample_1(p, L);

    j_re[0] = 1;
    for(int i = 1; i < L; i++){
        j_re[i] = 0;
    }
    
    uint64_t ret = 0;
    commit(b, r_e, r, c, m, d, r_f, r_g, C_upp, D_upp, E_upp, F_upp, G_0, W_upp, G_1, ConstTPrime, a, j_re, y, gamma, mat1, mat2, H_0_HAT, H_1_HAT, NHAT_mat);
    response(f, f_J, z_a, z, z_J, a, b, c, d, j_re, m, r, r_e, r_f, r_g, x);
    ret = verify(p, C_upp, D_upp, E_upp, F_upp, G_0, G_1, ConstTPrime, W_upp, x, y, gamma, f, f_J, z_a, z, z_J, mat1, mat2, H_0_HAT, H_1_HAT, NHAT_mat);
    
    //not test until i upload new code
    if(ret == 0){
        printf("Verified!!!");
    }

    return 0;

}
