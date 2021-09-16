#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "poly_param.h"
#include "poly_q.h"
#include "poly_red.h"
#include "comp.h"

//(398595837, 4294966337)
static const uint64_t omega_q[R] = {0, 1976180467LL, 973052952LL, -1284852282LL, 1945043797LL, 1291445467LL, -1777304198LL, -1712841977LL, -274282460LL, 901200135LL, 1480618753LL, -803283922LL, -1519077718LL, 1542850767LL, -942335666LL, -741384793LL, 1403690LL, 64886135LL, -451317068LL, 1714647223LL, -1794410942LL, 1786707388LL, -247775179LL, -1431557145LL, 1699289203LL, 1509901965LL, -1379590971LL, -1392633560LL, -999616345LL, 1163345716LL, 1978912114LL, -1029847677LL};

//static const uint64_t perm_ntt[R] = {0, 32, 8, 44, 5, 39, 13, 40, 2, 35, 10, 47, 7, 37, 15, 43, 1, 33, 9, 45, 4, 38, 12, 41, 3, 34, 11, 46, 6, 36, 14, 42, 63, 31, 55, 19, 58, 24, 50, 23, 61, 28, 53, 16, 56, 26, 48, 20, 62, 30, 54, 18, 59, 25, 51, 22, 60, 29, 52, 17, 57, 27, 49, 21};

//static const uint64_t perm_inv_aut[D] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198, 201, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231, 234, 237, 240, 243, 246, 249, 252, 255, 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 95, 98, 101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191, 194, 197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233, 236, 239, 242, 245, 248, 251, 254, 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 100, 103, 106, 109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154, 157, 160, 163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 202, 205, 208, 211, 214, 217, 220, 223, 226, 229, 232, 235, 238, 241, 244, 247, 250, 253};

static const uint64_t omega_qhat[R] = {0, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, 12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL, -12LL};

/* in bit-reversed order */
static inline void ntt_q_bitrev(POLY_Q *a)
{
	uint64_t t = D >> 1, m, i, j, j_first, j_last;
	uint64_t u, v;

	/* implicitly convert to Montgomery form and perform the butterfly operations */
	for (i = 0; i < D >> 1; i++)
	{
		u = montgomery_q(a->poly[i], MONTGOMERY_CONVERT_FACTOR_Q);
		v = montgomery_q(a->poly[i + (D >> 1)], omega_q[1]);
		
		a->poly[i] = con_sub(u + v, Q);
		a->poly[i + (D >> 1)] = con_add(u - v, Q);
	}

	for (m = 2; m < R; m <<= 1)
	{
		t >>= 1;
		for (i = 0; i < m; i++)
		{
			j_first = 2 * i * t;
			j_last = j_first + t;
			/* butterfly: (u,v)-->(u+v*omega,u-v*omega) */
			for (j = j_first; j < j_last; j++)
			{
				u = a->poly[j];
				v = montgomery_q(a->poly[j + t], omega_q[m + i]);
				
				a->poly[j] = con_sub(u + v, Q);
				a->poly[j + t] = con_add(u - v, Q);
			}
		}
	}
}

/* in the ordering (X - w^{(-1)^s * 3^r % 128}) for s = {0..1}, r = {0..31} */
void ntt_q(POLY_Q *a)
{
	POLY_Q tmp;
	
	uint64_t i, j;
	
	memcpy(&tmp, a, sizeof(POLY_Q));
	ntt_q_bitrev(&tmp);
	
	for (i = 0; i < R; i++)
	{
		for (j = 0; j < SUBRING_SIZE; j++)
		{
			a->poly[i * SUBRING_SIZE + j] = tmp.poly[i * SUBRING_SIZE + j];
		}
	}
}


/* in bit-reversed order (seems there is no need to permute) */
void ntt_qhat(POLY_QHAT *a)
{
	uint64_t t = D >> 1, m, i, j, j_first, j_last;
	uint64_t u, v;

	/* implicitly convert to Montgomery form and perform the butterfly operations */
	for (i = 0; i < D >> 1; i++)
	{
		u = montgomery_qhat(a->poly[i], MONTGOMERY_CONVERT_FACTOR_QHAT);
		v = montgomery_qhat(a->poly[i + (D >> 1)], omega_qhat[1]);
		
		a->poly[i] = con_sub(u + v, QHAT);
		a->poly[i + (D >> 1)] = con_add(u - v, QHAT);
	}

	for (m = 2; m < R; m <<= 1)
	{
		t >>= 1;
		for (i = 0; i < m; i++)
		{
			j_first = 2 * i * t;
			j_last = j_first + t;
			/* butterfly: (u,v)-->(u+v*omega,u-v*omega) */
			for (j = j_first; j < j_last; j++)
			{
				u = a->poly[j];
				v = montgomery_qhat(a->poly[j + t], omega_qhat[m + i]);
				
				a->poly[j] = con_sub(u + v, QHAT);
				a->poly[j + t] = con_add(u - v, QHAT);
			}
		}
	}
}


