#include <stdint.h>

#include "sammat.h"
#include "poly_q.h"
#include "fastrandombytes.h"
#include "param.h"
#include "littleendian.h"

#define SAMMAT_Q_BYTES 4LL
#define SAMMAT_Q_BOUND 4194256025LL

#define SAMMAT_QHAT_BYTES 1LL
// #define SAMMAT_QHAT_BOUND 4194256025LL

#define SAMMAT_MAT1_LEN 297216LL
#define SAMMAT_MAT2_LEN 2322LL
#define SAMMAT_MAT3_LEN 774LL


/* Mat_1 is n * (2 * l * \hat{n} * k), and Mat_2 is 1 * m */

/* sample Mat_1 matrix */
static inline void sample_mat_g(POLY_Q out[][N], const uint64_t m, const uint64_t rej_len)
{
	static unsigned char r[SAMMAT_MAT1_LEN * SAMMAT_Q_BYTES];
	uint64_t i, j, k, x;
	
	unsigned char *r_head = r;
	
	fastrandombytes_setiv_mat1();
	
	fastrandombytes_pub(r, SAMMAT_Q_BYTES * rej_len);
	
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < D; k++)
			{
				do
				{
					x = load_32(r_head);
					r_head += SAMMAT_Q_BYTES;
				} while (x >= SAMMAT_Q_BOUND); 
				
				out[i][j].poly[k] = x % Q;
			}
		}
	}
}

void sample_mat1(POLY_Q out[][N])
{
	sample_mat_g(out, 2 * L * NHAT * K, SAMMAT_MAT1_LEN);
}


/* sample Mat_2 matrix */
void sample_mat2(POLY_Q out[][1])
{
	static unsigned char r[SAMMAT_MAT2_LEN * SAMMAT_Q_BYTES];
	uint64_t i, j, k, x;
	
	unsigned char *r_head = r;
	
	fastrandombytes_setiv_mat2();
	
	fastrandombytes_pub(r, SAMMAT_Q_BYTES * SAMMAT_MAT2_LEN);
	
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < 1; j++)
		{
			for (k = 0; k < D; k++)
			{
				do
				{
					x = load_32(r_head);
					r_head += SAMMAT_Q_BYTES;
				} while (x >= SAMMAT_Q_BOUND); 
				
				out[i][j].poly[k] = x % Q;
			}
		}
	}
}

// /* --------------------
//  * |I_n | Mat_3       |
//  * |------------------|
//  * 
//  * where Mat_3 is n * (2 * l * \hat{n} * k) */
// /* sample Mat_3 matrix */
// void sample_mat3(POLY_Q out[][N])
// {
// 	static unsigned char r[SAMMAT_MAT3_LEN * SAMMAT_Q_BYTES];
// 	uint64_t i, j, k, x;
	
// 	unsigned char *r_head = r;
	
// 	fastrandombytes_setiv_mat3();
	
// 	fastrandombytes_pub(r, SAMMAT_Q_BYTES * SAMMAT_MAT3_LEN);
	
// 	for (i = 0; i < 2 * L * NHAT * K; i++)
// 	{
// 		for (j = 0; j < N; j++)
// 		{
// 			for (k = 0; k < D; k++)
// 			{
// 				do
// 				{
// 					x = load_32(r_head);
// 					r_head += SAMMAT_Q_BYTES;
// 				} while (x >= SAMMAT_Q_BOUND); 
				
// 				out[i][j].poly[k] = x % Q;
// 			}
// 		}
// 	}
// }


/* sample Mat_3 matrix which is \hat{n} * (\hat{n} * k) */
void sample_mat3(POLY_QHAT out[][NHAT])
{
	static unsigned char r[SAMMAT_MAT3_LEN * SAMMAT_QHAT_BYTES];
	uint64_t i, j, k, x;
	
	unsigned char *r_head = r;
	
	fastrandombytes_setiv_mat3();
	
	fastrandombytes_pub(r, SAMMAT_QHAT_BYTES * SAMMAT_MAT3_LEN);
	
	for (i = 0; i < NHAT * K; i++)
	{
		for (j = 0; j < NHAT; j++)
		{
			for (k = 0; k < D; k++)
			{
				x = load_8(r_head);
				r_head += SAMMAT_QHAT_BYTES;
				out[i][j].poly[k] = x % QHAT;
			}
		}
	}
}
