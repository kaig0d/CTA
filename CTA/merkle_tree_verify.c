#include<poly_q.h>
#include<poly_red.h>
#include<param.h>
#include<math.h>
#include<poly_mult.h>
#include<comp.h>

int verify(POLY_QHAT *p, POLY_Q *C_upp, POLY_Q *D_upp, POLY_Q *E_upp, POLY_Q *F_upp, POLY_Q *G_0, POLY_Q *G_1, POLY_Q *ConstTPrime, POLY_Q *W_upp, POLY_Q *x, POLY_Q *y, POLY_Q * gamma, POLY_Q *f,  POLY_Q *f_J, POLY_Q *z_a, POLY_Q *z, POLY_Q *z_J, POLY_Q mat1[][N], POLY_Q mat2[][1], POLY_Q H_0_HAT[][L * NHAT], POLY_Q H_1_HAT[][L * NHAT], POLY_Q NHAT_mat[][L * NHAT]){
    POLY_Q *f_ntt;
    POLY_Q *z_a_ntt;
    POLY_Q *z_ntt;
    POLY_Q *f_J_ntt;
    POLY_Q *z_J_ntt;

    POLY_Q *D_x_C;
    POLY_Q *A_f_B_z_a;
    POLY_Q *A_f;
    POLY_Q *B_z_a;
    POLY_Q *F_x_E;
    POLY_Q *A_1_f_J_B_1_z_J;
    POLY_Q *A_1_f_J;
    POLY_Q *B_1_z_J;
    POLY_Q *W_x_G_0;
    POLY_Q *K_0_z;
    POLY_Q *K_1_z;
    POLY_Q *x_power2;

    POLY_R fJ_x_fJ[L];
    POLY_Q *X_G1;
    POLY_Q *x_G1_K1_z;
    POLY_Q sum_1[L];
    POLY_Q *sum_1_plus;
    POLY_Q x_fJ[L];
    POLY_Q sum_2[NHAT * L];
    POLY_Q *sum_2_plus;
    POLY_Q sum_3[NHAT * L];
    POLY_Q *sum_3_plus;
    POLY_Q *sum_1_sum_2;

    POLY_Q H_0_HAT_mul[2 * L * NHAT * K][L * NHAT];
    POLY_Q H_1_HAT_mul[2 * L * NHAT * K][L * NHAT];
    POLY_Q N_HAT_mul[2 * L * NHAT * K][L * NHAT];
    POLY_Q Mat_sum[2 * L * NHAT * K][L * NHAT];
    POLY_Q Mat_sum_f[L * NHAT];
    uint64_t z_2norm;
    uint64_t z_a2norm;
    uint64_t z_J2norm;
    uint64_t f_J2norm;
    uint64_t f_2norm;
    uint64_t ret;
    
    // uint64_t B_z;
    // uint64_t B_J;
    // uint64_t B_f;

    // B_z = 1.32 * D * sqrt(3 * M * D) * D * sqrt(3 * M * D) * sqrt(2 * M * D);
    // B_z = B_z * B_z;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < D; j++){
            z_2norm += (int64_t)(z[i].poly[j]) * (int64_t)(z[i].poly[j]);
            z_a2norm += (int64_t)(z_a[i].poly[j]) * (int64_t)(z_a[i].poly[j]);
            z_J2norm += (int64_t)(z_J[i].poly[j]) * (int64_t)(z_J[i].poly[j]);
        }
    }

    if((z_2norm > B_z) || (z_a2norm > B_z) || (z_J2norm > B_z)){
        return 1;
    }

    // B_J = 1.32 * D * sqrt(L) * D * sqrt(L) * sqrt(2 * L * D);
    // B_J = B_J * B_J;
    //power and poly mul 
    for(int i = 0; i < L; i++){
        for(int j = 0; j < D; j++){
            f_J2norm += (int64_t)(f_J[i].poly[j]) * (int64_t)(f_J[i].poly[j]);

        }
    }

    if(f_J2norm > B_J){
        return 1;
    }

    // B_f = 0.675 * D * sqrt(2 * L * NHAT * K * D) * D * sqrt(2 * L * NHAT * K * D) * sqrt(L * NHAT * K * D);
    // B_f = B_f * B_f;

    for(int i = 0; i < 2 * L * NHAT * K; i++){
        for(int j = 0; j < D; j++){
            f_2norm += (int64_t)(f[i].poly[j]) * (int64_t)(f[i].poly[j]);

        }
    }

    if(f_2norm > B_f){
        return 1;
    }
    // D + x*C
    for (int i = 0; i < N; i++)
    {
        mult_r(D_x_C + i, x, C_upp + i);
        
        for (int j = 0; j < D; j++)
        {
            D_x_C[i].poly[j] = con_sub(D_x_C[i].poly[j] + D_upp[i].poly[j], Q);
        }
    }

    /* ntt */
	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			f_ntt[i].poly[j] = con_add(f[i].poly[j], Q);
		}
		
		ntt_q(f_ntt + i);
	}
    /*Af = A * f*/
    for (int i = 0; i < N; i++)
	{
		memcpy(A_f + i, f_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < 2 * L * NHAT * K - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(A_f + j, mat1[i] + j, f_ntt + N + i);
		}
	}

    /* ntt */
	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			z_a_ntt[i].poly[j] = con_add(z_a[i].poly[j], Q);
		}
		
		ntt_q(z_a_ntt + i);
	}
    /* Bza= B * za*/
    for (int i = 0; i < N; i++)
	{
		memcpy(B_z_a + i, z_a_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(B_z_a + j, mat1[i] + j, z_a_ntt + N + i);
		}
	}

    // A_f_B_z_a = A_f + B_z_a
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            A_f_B_z_a[i].poly[j] = con_sub(A_f[i].poly[j] + B_z_a[i].poly[j], Q);

        }
    }

    ret = 0; 
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            ret += ct_eq(A_f_B_z_a[i].poly[j], D_x_C[i].poly[j]);
        }

    }

    if(ret != D){
        return 1;
    }

    // F + x*E
    for (int i = 0; i < N; i++)
    {
        mult_r(F_x_E + i, x, E_upp + i);
        
        for (int j = 0; j < D; j++)
        {
            F_x_E[i].poly[j] = con_sub(F_x_E[i].poly[j] + F_upp[i].poly[j], Q);
        }
    }

    //
    /* ntt */
	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < D; j++)
		{
			f_J_ntt[i].poly[j] = con_add(f_J[i].poly[j], Q);
		}
		
		ntt_q(f_J_ntt + i);
	}
    /*A_1_f_J = A1 * fJ*/
    for (int i = 0; i < N; i++)
	{
		memcpy(A_1_f_J + i, f_J_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < L - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(A_1_f_J + j, mat1[i] + j, f_J_ntt + N + i);
		}
	}

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			z_J_ntt[i].poly[j] = con_add(z_J[i].poly[j], Q);
		}
		
		ntt_q(z_J_ntt + i);
	}
    /* B_1_z_J= B1 * zJ*/
    for (int i = 0; i < N; i++)
	{
		memcpy(B_1_z_J + i, z_J_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(B_1_z_J + j, mat1[i] + j, z_J_ntt + N + i);
		}
	}

    // A_1_f_J_B_1_z_J = A_1_f_J + B_1_z_J
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            A_1_f_J_B_1_z_J[i].poly[j] = con_sub(A_1_f_J[i].poly[j] + B_1_z_J[i].poly[j], Q);

        }
    }

    ret = 0; 
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            ret += ct_eq(A_1_f_J_B_1_z_J[i].poly[j], F_x_E[i].poly[j]);

        }

    }

    if(ret != D){
        return 1;
    }

    // W + x*G0
    for (int i = 0; i < N; i++)
    {
        mult_r(W_x_G_0 + i, x, G_0 + i);
        
        for (int j = 0; j < D; j++)
        {
            W_x_G_0[i].poly[j] = con_sub(W_x_G_0[i].poly[j] + W_upp[i].poly[j], Q);
        }
    }

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			z_ntt[i].poly[j] = con_add(z[i].poly[j], Q);
		}
		
		ntt_q(z_ntt + i);
	}
    /*K_0_z = K0 * z*/
    for (int i = 0; i < N; i++)
	{
		memcpy(K_0_z + i, z_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(K_0_z + j, mat1[i] + j, z_ntt + N + i);
		}
	}


    ret = 0; 
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            ret += ct_eq(W_x_G_0[i].poly[j], K_0_z[i].poly[j]);

        }
  
    }

    if(ret != D){
        return 1;
    }


    // \sum_{i=1}^{l} ri(x - fJ[i])fJ[i]
    for(int i = 0; i < L; i++){
        mult_fxf(fJ_x_fJ + i, f_J + i, x);
        mult_r(sum_1 + i, gamma + i, fJ_x_fJ + i);
    }

    for(int i = 0; i < L; i++){
        for(int j = 0; j < D; j++){
            sum_1_plus->poly[j] = con_sub(sum_1_plus->poly[j] + sum_1[i].poly[j], Q);
        }
    }
    
    //sum2
    for(int i = 0 ; i < L; i++){
        for(int j = i * 2 * NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                mult_r(H_0_HAT_mul[j]+k, H_0_HAT[j]+k, f_J + i);

            }
        }
    }

    for(int i = 0; i < L; i++){
        for(int j = 0; j < D; j++){
            x_fJ[i].poly[j] = con_add(x_fJ[i].poly[j] - f_J[i].poly[j], Q);
        }
    }

    for(int i = 0 ; i < L; i++){
        for(int j = i * 2 * NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                mult_r(H_1_HAT_mul[j]+k, H_1_HAT[j]+k, x_fJ + i);

            }
        }
    }

    //change
    for(int i = 0; i < 2 * L * NHAT * K; i++){
        for(int j = 0; j < L * NHAT; j++){
            mult_r(N_HAT_mul[i] + j, NHAT_mat[i] + j, x);

        }
    }

    for(int i = 0; i < 2 * L * NHAT * K; i++){
        for(int j = 0; j < L * NHAT; j++){
            for(int k = 0; k < D; k++){
                Mat_sum[i][j].poly[k] = con_sub(H_0_HAT_mul[i][j].poly[k] + H_1_HAT_mul[i][j].poly[k] - N_HAT_mul[i][j].poly[k], Q);
                Mat_sum[i][j].poly[k] = con_add(H_0_HAT_mul[i][j].poly[k] + H_1_HAT_mul[i][j].poly[k] - N_HAT_mul[i][j].poly[k], Q);
            }

        }
    }

	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < L * NHAT; j++)
		{
			mult_r(Mat_sum_f + j, Mat_sum[i] + j, f_ntt + i);
		}
	}

    for(int i = 0; i < L * NHAT; i++){
        mult_r(sum_2 + i, Mat_sum_f + i, y + i);

    }

    for(int i = 0; i < L * NHAT; i++){
        for(int j = 0; j < D; j++){
            sum_2_plus->poly[j] = con_sub(sum_2_plus->poly[j] + sum_2[i].poly[j], Q);
        }
    }

    for(int i = 0; i < D; i++){
        sum_1_sum_2->poly[i] = con_sub(sum_1_plus->poly[i] + sum_2_plus->poly[i], Q);

    }


    //x^{2}\sum_{i=0}^{\hat{n}*l-1}y[i]p[i]+(xG1-K1z)+ConsTPrime
    mult_rq(x_power2, x, x);

    //p is under R_{qhat}
    for(int i = 0; i < NHAT * L; i++){
        mult_rqhat(sum_3 + i, y + i, p +i);
        mult_r(sum_3 + i, sum_3 + i, x_power2);
    }

    for(int i = 0; i < NHAT * L; i++){
        for(int j = 0; j < D; j++){
            sum_3_plus->poly[j] = con_sub(sum_3_plus->poly[j] + sum_3[i].poly[j], Q);
        }
    }
    
    for(int i = 0; i < 1; i++){
        mult_r(X_G1 + i, x, G_1 + i);
    }

    // K1 * z
    for (int i = 0; i < 1; i++)
	{
		memcpy(K_1_z + i, z_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - 1; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			mult_plus_rq(K_1_z + j, mat2[i] + j, z_ntt + 1 + i);
		}
	}

    for(int i = 0; i < D; i++){
        x_G1_K1_z->poly[i] = con_add(X_G1->poly[i] - K_1_z->poly[i], Q);

    }

    for(int i = 0; i < D; i++){
        sum_3_plus->poly[i] = con_sub(sum_3_plus->poly[i] + x_G1_K1_z->poly[i] + ConstTPrime->poly[i], Q);
        sum_3_plus->poly[i] = con_add(sum_3_plus->poly[i] + x_G1_K1_z->poly[i] + ConstTPrime->poly[i], Q);
    }

    ret = 0; 
    for(int i = 0; i < D; i++){
        ret += ct_eq(sum_1_sum_2->poly[i], sum_3_plus->poly[i]);

    }
    if(ret != D){
        return 1;
    }


    return 0;

}