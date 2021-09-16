#include<param.h>
#include<poly_q.h>
#include<poly_red.h>
#include<sample.h>
#include<gaussian_avx.h>
#include<poly_mult.h>
#include<comp.h>

void commit(POLY_Q *b, POLY_Q *r_e, POLY_Q *r, POLY_Q *c, POLY_Q *m, POLY_Q *d, POLY_Q *r_f, POLY_Q *r_g, POLY_Q *C_upp, POLY_Q *D_upp, POLY_Q *E_upp, POLY_Q *F_upp, POLY_Q *G_0, POLY_Q *W_upp, POLY_Q *G_1, POLY_Q *ConstTPrime, POLY_Q *a, uint64_t j_re[L], POLY_Q  *y, POLY_Q *gamma, POLY_Q mat1[][N], POLY_Q mat2[][1], POLY_Q H_0_HAT[][L * NHAT], POLY_Q H_1_HAT[][L * NHAT], POLY_Q NHAT_mat[][L * NHAT]){
    POLY_Q *b_ntt;
	POLY_Q *r_e_ntt;
	POLY_Q *r_ntt;
	POLY_Q *c_ntt;
	POLY_Q *m_ntt;
	POLY_Q *d_ntt;
	POLY_Q *r_f_ntt;
	POLY_Q *r_g_ntt;
	POLY_Q *a_ntt;

	POLY_Q *C_upp_temp;
    POLY_Q *D_upp_temp;
    POLY_Q *E_upp_temp;
    POLY_Q *F_upp_temp;
    POLY_Q *PrimTPrime;
	POLY_Q *PrimT;
	POLY_Q *ConsT;
	POLY_Q j_H_0_HAT[2 * NHAT * L * K][L * NHAT];
	POLY_Q jBAR_H_1_HAT[2 * NHAT * L * K][L * NHAT];
	POLY_Q m_H_0_HAT[2 * NHAT * L * K][L * NHAT];
	POLY_Q m_H_1_HAT[2 * NHAT * L * K][L * NHAT];
	POLY_Q Mat_sum[2 * NHAT * L * K][L * NHAT];
	POLY_Q Mat_sum_c[L * NHAT];
	POLY_Q Mat_minus[2 * NHAT * L * K][L * NHAT];
	POLY_Q Mat_minus_a[L * NHAT];
	POLY_Q Mat_minus_c[L * NHAT];
	uint64_t j_re_BAR[L];
	POLY_Q gamma_m[L];
	POLY_Q *gamma_m_sum;
	POLY_Q gamma_m_power[L];
	POLY_Q *gamma_m_power_sum;
	POLY_Q K_1_r_g[1];
    

    sample_1(b, M);
    sample_1(r_e, M);
    sample_1(r, M);

    sample_b1(c);

    sample_b2(m);

    sample_b3(d);
    sample_b3(r_f);
    sample_b3(r_g);
    
    
    /* ntt */
	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			a_ntt[i].poly[j] = con_add(a[i].poly[j], Q);
		}
		
		ntt_q(a_ntt + i);
	}
    /*C = A * a*/
    for (int i = 0; i < N; i++)
	{
		memcpy(C_upp + i, a_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < 2 * L * NHAT * K - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(C_upp + j, mat1[i] + j, a_ntt + N + i);
		}
	}

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			b_ntt[i].poly[j] = con_add(b[i].poly[j], Q);
		}
		
		ntt_q(b_ntt + i);
	}

    /*C_temp = B * b*/
    for (int i = 0; i < N; i++)
	{
		memcpy(C_upp_temp + i, b_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(C_upp_temp + j, mat1[i] + j, b_ntt + N + i);
		}
	}

    /*C = C + C_temp*/
    for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			C_upp[i].poly[j] = con_sub(C_upp[i].poly[j] + C_upp_temp[i].poly[j], Q);
		}
	}

    /* ntt */
	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < D; j++)
		{
			c_ntt[i].poly[j] = con_add(c[i].poly[j], Q);
		}
		
		ntt_q(c_ntt + i);
	}

    /*D = A * c*/
    for (int i = 0; i < N; i++)
	{
		memcpy(D_upp + i, c_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < 2 * L * NHAT * K - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(D_upp + j, mat1[i] + j, c_ntt + N + i);
		}
	}

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			d_ntt[i].poly[j] = con_add(d[i].poly[j], Q);
		}
		
		ntt_q(d_ntt + i);
	}

    /*D_temp = B * d*/
    for (int i = 0; i < N; i++)
	{
		memcpy(D_upp_temp + i, d_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(D_upp_temp + j, mat1[i] + j, d_ntt + N + i);
		}
	}

    /*D = D + D_temp*/
    for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			D_upp[i].poly[j] = con_sub(D_upp[i].poly[j] + D_upp_temp[i].poly[j], Q);
		}
	}
	// j[i] is 0 or 1 
    /*E = A1 * j*/
	for (int i = 0; i < L; i++)
	{
		for(int j = 0; j < N; j++){
			if(j_re[i] == 0){
				for(int k = 0; k < D; k++){
					E_upp[j].poly[k] = 0;
				}
			}else{
				for(int m_indx = 0; m_indx < L; m_indx++){
					for(int k = 0; k < D; k++){
						//something wrong
						E_upp[j].poly[k] = con_sub(E_upp[j].poly[k] + mat1[m_indx][j].poly[k], Q);
					}
				}
			}
		}
	}

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			r_e_ntt[i].poly[j] = con_add(r_e[i].poly[j], Q);
		}
		
		ntt_q(r_e_ntt + i);
	}

    /*E_temp = B1 * r_e*/
    for (int i = 0; i < N; i++)
	{
		memcpy(E_upp_temp + i, r_e_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(E_upp_temp + j, mat1[i] + j, r_e_ntt + N + i);
		}
	}

    /*E = E + E_temp*/
    for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			E_upp[i].poly[j] = con_sub(E_upp[i].poly[j] + E_upp_temp[i].poly[j], Q);
		}
	}
    
    /* ntt */
	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < D; j++)
		{
			m_ntt[i].poly[j] = con_add(m[i].poly[j], Q);
		}
		
		ntt_q(m_ntt + i);   
	}
    /*F = A1 * m*/
    for (int i = 0; i < N; i++)
	{
		memcpy(F_upp + i, m_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < L - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(F_upp + j, mat1[i] + j, m_ntt + N + i);
		}
	}

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			r_f_ntt[i].poly[j] = con_add(r_f[i].poly[j], Q);
		}
		
		ntt_q(r_f_ntt + i);
	}

    /*F_temp = B1 * r_f*/
    for (int i = 0; i < N; i++)
	{
		memcpy(F_upp_temp + i, r_f_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(F_upp_temp + j, mat1[i] + j, r_f_ntt + N + i);
		}
	}

    /*F = F + F_temp*/
    for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < D; j++)
		{
			F_upp[i].poly[j] = con_sub(F_upp[i].poly[j] + F_upp_temp[i].poly[j], Q);
		}
	}

    /* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			r_ntt[i].poly[j] = con_add(r[i].poly[j], Q);
		}
		
		ntt_q(r_ntt + i);
	}

    /*G0 = K0 * r*/
    for (int i = 0; i < N; i++)
	{
		memcpy(G_0 + i, r_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(G_0 + j, mat1[i] + j, r_ntt + N + i);
		}
	}
    
	/* ntt */
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < D; j++)
		{
			r_g_ntt[i].poly[j] = con_add(r_g[i].poly[j], Q);
		}
		
		ntt_q(r_g_ntt + i);
	}

    /*W = K0 * r_g*/
    for (int i = 0; i < N; i++)
	{
		memcpy(W_upp + i, r_g_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mult_plus_rq(W_upp + j, mat1[i] + j, r_g_ntt + N + i);
		}
	}

	//PrimTPrime,ConstPrime

	//j is reversed
	//j_H_0_HAT
	for(int i = 0 ; i < L; i++){
		if(j_re[i] == 0){
			for(int j = i * 2 * NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
				for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
					for(int m_indx = 0; m_indx < D; m_indx++){
						j_H_0_HAT[j][k].poly[m_indx] = 0;
					}
					
				}
        	}
		}
        
    }

	for(int i = 0; i < L; i++){
		if(j_re[i] = 0){
			j_re_BAR[i] = 1;

		}else{
			j_re_BAR[i] = 0;
		}
	}

	//jBAR_H_1_HAT
	for(int i = 0 ; i < L; i++){
		if(j_re_BAR[i] == 0){
			for(int j = i * 2 * NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
				for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
					for(int m_indx = 0; m_indx < D; m_indx++){
						jBAR_H_1_HAT[j][k].poly[m_indx] = 0;
					}
					
				}
       		}
		}
        
    }

	//m_H_0_HAT
	for(int i = 0 ; i < L; i++){
        for(int j = i * 2 * NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                mult_r(m_H_0_HAT[j]+k, H_0_HAT[j]+k, m + i);

            }
        }
    }
	//m_H_1_HAT
	for(int i = 0 ; i < L; i++){
        for(int j = i * 2 * NHAT * K; j < (i + 1) * 2 * NHAT * K; j++){
            for(int k = i * NHAT; k < (i + 1) * NHAT; k++){
                mult_r(m_H_1_HAT[j]+k, H_1_HAT[j]+k, m + i);

            }
        }
    }

	for(int i = 0; i < 2 * L * NHAT * K; i++){
        for(int j = 0; j < L * NHAT; i++){
            for(int k = 0; k < D; k++){
                Mat_sum[i][j].poly[k] = con_sub(j_H_0_HAT[i][j].poly[k] + jBAR_H_1_HAT[i][j].poly[k] - NHAT_mat[i][j].poly[k], Q);
				Mat_sum[i][j].poly[k] = con_add(j_H_0_HAT[i][j].poly[k] + jBAR_H_1_HAT[i][j].poly[k] - NHAT_mat[i][j].poly[k], Q);
            }
        }
    }

	for(int i = 0; i < 2 * L * NHAT * K; i++){
        for(int j = 0; j < L * NHAT; i++){
            for(int k = 0; k < D; k++){
                Mat_minus[i][j].poly[k] = con_add(m_H_0_HAT[i][j].poly[k] - m_H_1_HAT[i][j].poly[k], Q);

            }

        }
    }

	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < L * NHAT; j++)
		{
			mult_plus_rq(Mat_sum_c + j, Mat_sum[i] + j, c_ntt + i);
		}
	}


	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < L * NHAT; j++)
		{
			//something wrong
			mult_plus_rq(Mat_minus_a + j, Mat_minus[i] + j, a_ntt + i);
		}
	}

	//PrimT
	for(int i = 0; i < L * NHAT; i++){
		mult_plus_rq(PrimT, Mat_sum_c + i, y + i);
	}

	for(int i = 0; i < L * NHAT; i++){
		for(int j = 0; j < D; j++){
			PrimT->poly[j] = con_sub(PrimT->poly[j] + Mat_minus_a->poly[j], Q);
		}
	}

	for (int i = 0; i < 2 * L * NHAT * K; i++)
	{
		for (int j = 0; j < L * NHAT; j++)
		{
			mult_plus_rq(Mat_minus_c + j, Mat_minus[i] + j, c_ntt + i);
		}
	}

	//ConsT, NTT domain?
	for(int i = 0; i < L * NHAT; i++){
		mult_plus_rq(ConsT, Mat_minus_c + i, y + i);
	}


	for(int i = 0; i < L; i++){
		for(int j = 0; j < D; j++){
			gamma_m[i].poly[j] = con_sub(gamma[i].poly[j] * m[i].poly[j], Q);
		}
	}

	for(int i = 0; i < L; i++){
		for(int j = 0; j < D; j++){
			gamma_m_sum->poly[j] = con_sub(gamma_m_sum->poly[j] + gamma_m[i].poly[j] *(1 - 2 * j_re[i]), Q); 
		}
		
	}
	//PrimTPrime
	for(int i = 0; i < D; i++){
		PrimTPrime->poly[i] = con_sub(PrimT->poly[i] + gamma_m_sum->poly[i], Q);
	}


    /*G1 = K1 * r*/
    for (int i = 0; i < 1; i++)
	{
		memcpy(G_1 + i, r_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - 1; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			mult_plus_rq(G_1 + j, mat2[i] + j, r_ntt + 1 + i);
		}
	}

    /*G1 = G1 + PrimTPrime*/
    for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < D; j++)
		{
			G_1[i].poly[j] = con_sub(G_1[i].poly[j] + PrimTPrime[i].poly[j], Q);
		}
	}

	for(int i = 0; i < L; i++){
		for(int j = 0; j < D; j++){
			gamma_m_power[i].poly[j] = con_sub(gamma_m[i].poly[j] * m[i].poly[j], Q);
		}
	}

	for(int i = 0; i < L; i++){
		for(int j = 0; j < D; j++){
			gamma_m_power_sum->poly[j] = con_sub(gamma_m_power_sum->poly[j] + gamma_m_power[i].poly[j], Q); 
		}
		
	}

    /*G1 = K1 * r*/
    for (int i = 0; i < 1; i++)
	{
		memcpy(K_1_r_g + i, r_g_ntt + i, sizeof(POLY_Q));
	}

	for (int i = 0; i < M - 1; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			mult_plus_rq(K_1_r_g + j, mat2[i] + j, r_g_ntt + 1 + i);
		}
	}

	for(int i = 0; i < D; i++){
		ConstTPrime->poly[i] = con_sub(ConsT->poly[i] - gamma_m_power_sum->poly[i] + K_1_r_g->poly[i], Q);
		ConstTPrime->poly[i] = con_add(ConsT->poly[i] - gamma_m_power_sum->poly[i] + K_1_r_g->poly[i], Q);
	}
	
}