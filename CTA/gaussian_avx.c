#include "gaussian_avx.h"
#include <stdint.h>
#include "fastrandombytes.h"
#include <x86intrin.h>
#include "randombytes.h"
#include "param.h"
#include "littleendian.h"
#include "poly_q.h"
#include "comp.h"
#include <math.h>

#define CDT_ENTRY_SIZE 16
#define CDT_LOW_MASK 0x7fffffffffffffff
#define CDT_LENGTH 9 /* [0..tau*sigma]=[0..9] */

#define BERNOULLI_ENTRY_SIZE 9 /* 72bit randomness */

/* -1/k^2 */
#define BINARY_SAMPLER_B1_K_2_INV (-1.0/(B1_K * B1_K))
#define BINARY_SAMPLER_B2_K_2_INV (-1.0/(B2_K * B2_K))
#define BINARY_SAMPLER_B3_K_2_INV (-1.0/(B3_K * B3_K))

#define EXP_MANTISSA_PRECISION 52
#define EXP_MANTISSA_MASK ((1LL << EXP_MANTISSA_PRECISION) - 1)
#define R_MANTISSA_PRECISION (EXP_MANTISSA_PRECISION + 1)
#define R_MANTISSA_MASK ((1LL << R_MANTISSA_PRECISION) - 1)
#define R_EXPONENT_L (8 * BERNOULLI_ENTRY_SIZE - R_MANTISSA_PRECISION)

#define DOUBLE_ONE (1023LL << 52)

#define UNIFORM_SIZE_B1_K 2
#define UNIFORM_REJ_B1_K 20
#define BARRETT_BITSHIFT_B1_K (UNIFORM_SIZE_B_K * 8)
#define BARRETT_FACTOR_B1_K ((1LL << BARRETT_BITSHIFT_B1_K) / B1_K)
#define UNIFORM_Q_B1_K (B1_K * BARRETT_FACTOR_B1_K)

#define UNIFORM_SIZE_B2_K 2
#define UNIFORM_REJ_B2_K 20
#define BARRETT_BITSHIFT_B2_K (UNIFORM_SIZE_B2_K * 8)
#define BARRETT_FACTOR_B2_K ((1LL << BARRETT_BITSHIFT_B2_K) / B2_K)
#define UNIFORM_Q_B2_K (B2_K * BARRETT_FACTOR_B2_K)

#define UNIFORM_SIZE_B3_K 3
#define UNIFORM_REJ_B3_K 30
#define BARRETT_BITSHIFT_B3_K (UNIFORM_SIZE_B3_K * 8)
#define BARRETT_FACTOR_B3_K ((1LL << BARRETT_BITSHIFT_B3_K) / B3_K)
#define UNIFORM_Q_B3_K (B3_K * BARRETT_FACTOR_B3_K)

#define UNIFORM_R_SIZE_MAX 90

#define BASE_TABLE_SIZE (4 * CDT_ENTRY_SIZE)
#define BERNOULLI_TABLE_SIZE (4 * BERNOULLI_ENTRY_SIZE)

/* CDT table */
static const __m256i V_CDT[][2] = {{{2200310400551559144, 2200310400551559144, 2200310400551559144, 2200310400551559144}, {3327841033070651387, 3327841033070651387, 3327841033070651387, 3327841033070651387}},
{{7912151619254726620, 7912151619254726620, 7912151619254726620, 7912151619254726620}, {380075531178589176, 380075531178589176, 380075531178589176, 380075531178589176}},
{{5167367257772081627, 5167367257772081627, 5167367257772081627, 5167367257772081627}, {11604843442081400, 11604843442081400, 11604843442081400, 11604843442081400}},
{{5081592746475748971, 5081592746475748971, 5081592746475748971, 5081592746475748971}, {90134450315532, 90134450315532, 90134450315532, 90134450315532}},
{{6522074513864805092, 6522074513864805092, 6522074513864805092, 6522074513864805092}, {175786317361, 175786317361, 175786317361, 175786317361}},
{{2579734681240182346, 2579734681240182346, 2579734681240182346, 2579734681240182346}, {85801740, 85801740, 85801740, 85801740}},
{{8175784047440310133, 8175784047440310133, 8175784047440310133, 8175784047440310133}, {10472, 10472, 10472, 10472}},
{{2947787991558061753, 2947787991558061753, 2947787991558061753, 2947787991558061753}, {0, 0, 0, 0}},
{{22489665999543, 22489665999543, 22489665999543, 22489665999543}, {0, 0, 0, 0}}};

static const __m256i V_CDT_LOW_MASK = {CDT_LOW_MASK, CDT_LOW_MASK, CDT_LOW_MASK, CDT_LOW_MASK};

static const __m256i V_B1_K = {B1_K, B1_K, B1_K, B1_K};
static const __m256i V_B2_K = {B2_K, B2_K, B2_K, B2_K};
static const __m256i V_B3_K = {B3_K, B3_K, B3_K, B3_K};

/* coefficients of the exp evaluation polynomial */
static const __m256i EXP_COFF[] = {{0x3e833b70ffa2c5d4, 0x3e833b70ffa2c5d4, 0x3e833b70ffa2c5d4, 0x3e833b70ffa2c5d4},
								   {0x3eb4a480fda7e6e1, 0x3eb4a480fda7e6e1, 0x3eb4a480fda7e6e1, 0x3eb4a480fda7e6e1},
								   {0x3ef01b254493363f, 0x3ef01b254493363f, 0x3ef01b254493363f, 0x3ef01b254493363f},
								   {0x3f242e0e0aa273cc, 0x3f242e0e0aa273cc, 0x3f242e0e0aa273cc, 0x3f242e0e0aa273cc},
								   {0x3f55d8a2334ed31b, 0x3f55d8a2334ed31b, 0x3f55d8a2334ed31b, 0x3f55d8a2334ed31b},
								   {0x3f83b2aa56db0f1a, 0x3f83b2aa56db0f1a, 0x3f83b2aa56db0f1a, 0x3f83b2aa56db0f1a},
								   {0x3fac6b08e11fc57e, 0x3fac6b08e11fc57e, 0x3fac6b08e11fc57e, 0x3fac6b08e11fc57e},
								   {0x3fcebfbdff556072, 0x3fcebfbdff556072, 0x3fcebfbdff556072, 0x3fcebfbdff556072},
								   {0x3fe62e42fefa7fe6, 0x3fe62e42fefa7fe6, 0x3fe62e42fefa7fe6, 0x3fe62e42fefa7fe6},
								   {0x3ff0000000000000, 0x3ff0000000000000, 0x3ff0000000000000, 0x3ff0000000000000}};
								   
static const __m256d V_INT64_DOUBLE = {0x0010000000000000, 0x0010000000000000, 0x0010000000000000, 0x0010000000000000};
static const __m256d V_DOUBLE_INT64 = {0x0018000000000000, 0x0018000000000000, 0x0018000000000000, 0x0018000000000000};

static const __m256i V_EXP_MANTISSA_MASK = {EXP_MANTISSA_MASK, EXP_MANTISSA_MASK, EXP_MANTISSA_MASK, EXP_MANTISSA_MASK};
static const __m256i V_RES_MANTISSA = {1LL << EXP_MANTISSA_PRECISION, 1LL << EXP_MANTISSA_PRECISION, 1LL << EXP_MANTISSA_PRECISION, 1LL << EXP_MANTISSA_PRECISION};
static const __m256i V_RES_EXPONENT = {R_EXPONENT_L - 1023 + 1, R_EXPONENT_L - 1023 + 1, R_EXPONENT_L - 1023 + 1, R_EXPONENT_L - 1023 + 1};
static const __m256i V_R_MANTISSA_MASK = {R_MANTISSA_MASK, R_MANTISSA_MASK, R_MANTISSA_MASK, R_MANTISSA_MASK};
static const __m256i V_1 = {1, 1, 1, 1};
static const __m256i V_DOUBLE_ONE = {DOUBLE_ONE, DOUBLE_ONE, DOUBLE_ONE, DOUBLE_ONE};

static const __m256d V_B1_K_2_INV = {BINARY_SAMPLER_B1_K_2_INV, BINARY_SAMPLER_B1_K_2_INV, BINARY_SAMPLER_B1_K_2_INV, BINARY_SAMPLER_B1_K_2_INV};
static const __m256d V_B2_K_2_INV = {BINARY_SAMPLER_B2_K_2_INV, BINARY_SAMPLER_B2_K_2_INV, BINARY_SAMPLER_B2_K_2_INV, BINARY_SAMPLER_B2_K_2_INV};
static const __m256d V_B3_K_2_INV = {BINARY_SAMPLER_B3_K_2_INV, BINARY_SAMPLER_B3_K_2_INV, BINARY_SAMPLER_B3_K_2_INV, BINARY_SAMPLER_B3_K_2_INV};

static const __m128d REJ_V_B1_K_2_INV = {-BINARY_SAMPLER_B1_K_2_INV, 0};
static const __m128d REJ_V_B2_K_2_INV = {-BINARY_SAMPLER_B2_K_2_INV, 0};
static const __m128d REJ_V_B3_K_2_INV = {-BINARY_SAMPLER_B3_K_2_INV, 0};

#define REJ_BERNOULLI_ENTRY_SIZE 9 /* 64bit randomness */
#define REJ_R_EXPONENT_L (8 * REJ_BERNOULLI_ENTRY_SIZE - R_MANTISSA_PRECISION)
//how to calc
#define PHI_B1 0.675
#define PHI_B2 1.32
#define PHI_B3 1.32
#define MU_B1 (1.0 / (2.0 * PHI_B1 * PHI_B1) / M_LN2)
#define MU_B2 (1.0 / (2.0 * PHI_B2 * PHI_B2) / M_LN2)
#define MU_B3 (1.0 / (2.0 * PHI_B3 * PHI_B3) / M_LN2)

static const __m128d V_MU_B1 = {MU_B1, 0};
static const __m128d V_MU_B2 = {MU_B2, 0};
static const __m128d V_MU_B3 = {MU_B3, 0};

/* constant time CDT sampler */
static inline __m256i cdt_sampler(const unsigned char *r)
{
	__m256i x = _mm256_setzero_si256();
	__m256i r1, r2;
	__m256i r1_lt_cdt0, r2_lt_cdt1;
	__m256i r2_eq_cdt1;
	__m256i b;
	
	uint32_t i;
	
	r1 = _mm256_loadu_si256((__m256i *)r);
	r2 = _mm256_loadu_si256((__m256i *)(r + 32));
	
	r1 = _mm256_and_si256(r1, V_CDT_LOW_MASK);
	r2 = _mm256_and_si256(r2, V_CDT_LOW_MASK);

	for (i = 0; i < CDT_LENGTH; i++)
	{
		r1_lt_cdt0 = _mm256_sub_epi64(r1, V_CDT[i][0]);

		r2_lt_cdt1 = _mm256_sub_epi64(r2, V_CDT[i][1]);
		r2_eq_cdt1 = _mm256_cmpeq_epi64(r2, V_CDT[i][1]);

		b = _mm256_and_si256(r1_lt_cdt0, r2_eq_cdt1);
		b = _mm256_or_si256(b, r2_lt_cdt1);
		b = _mm256_srli_epi64(b, 63);

		x = _mm256_add_epi64(x, b);
	}

	return x;
}

/* constant time Bernoulli sampler
 * we directly compute exp(-x/(2*sigma^2)), 
 * since sigma_0=sqrt(1/2ln2), exp(-x/(2*sigma^2))=2^(-x/k^2), 
 * we use a polynomial to directly evaluate 2^(-x/k^2) */
static inline void bernoulli_sampler(uint64_t *b, __m256i x, const unsigned char *r, const __m256d v_k_2_inv)
{	
	__m256d vx, vx_1, vx_2, vsum;
	__m256i vt, k, vres, vres_mantissa, vres_exponent, vr_mantissa, vr_exponent, vr_exponent2, vres_eq_1, vr_lt_vres_mantissa, vr_lt_vres_exponent;

	/* 2^x=2^(floor(x)+a)=2^(floor(x))*2^a, where a is in [0,1]
	 * we only evaluate 2^a by using a polynomial */
	x = _mm256_or_si256(x, _mm256_castpd_si256(V_INT64_DOUBLE));
	vx = _mm256_sub_pd(_mm256_castsi256_pd(x), V_INT64_DOUBLE);
	vx = _mm256_mul_pd(vx, v_k_2_inv);
	
	vx_1 = _mm256_floor_pd(vx);
	vx_2 = _mm256_add_pd(vx_1, V_DOUBLE_INT64);
	vt = _mm256_sub_epi64(_mm256_castpd_si256(vx_2), _mm256_castpd_si256(V_DOUBLE_INT64));	
	vt = _mm256_slli_epi64(vt, 52);
	
	/* evaluate 2^a */
	vx_2 = _mm256_sub_pd(vx, vx_1);
	vsum = _mm256_fmadd_pd(_mm256_castsi256_pd(EXP_COFF[0]), vx_2, _mm256_castsi256_pd(EXP_COFF[1]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[2]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[3]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[4]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[5]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[6]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[7]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[8]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[9]));
	
	/* combine to compute 2^x */
	vres = _mm256_add_epi64(vt, _mm256_castpd_si256(vsum));

	/* compute the Bernoulli value */
	vres_mantissa = _mm256_and_si256(vres, V_EXP_MANTISSA_MASK);
	vres_mantissa = _mm256_or_si256(vres_mantissa, V_RES_MANTISSA);
	
	vres_exponent = _mm256_srli_epi64(vres, EXP_MANTISSA_PRECISION);
	vres_exponent = _mm256_add_epi64(vres_exponent, V_RES_EXPONENT);
	vres_exponent = _mm256_sllv_epi64(V_1, vres_exponent);
	
	vr_mantissa = _mm256_loadu_si256((__m256i *)r);
	vr_exponent = _mm256_srli_epi64(vr_mantissa, R_MANTISSA_PRECISION);
	vr_mantissa = _mm256_and_si256(vr_mantissa, V_R_MANTISSA_MASK);
	vr_exponent2 = _mm256_set_epi64x(r[35], r[34], r[33], r[32]);
	vr_exponent2 = _mm256_slli_epi64(vr_exponent2, 64 - R_MANTISSA_PRECISION);
	vr_exponent = _mm256_or_si256(vr_exponent, vr_exponent2);

	/* (res == 1.0) || ((r_mantissa < res_mantissa) && (r_exponent < (1 << res_exponent))) */
	vres_eq_1 = _mm256_cmpeq_epi64(vres, V_DOUBLE_ONE);
	vr_lt_vres_mantissa = _mm256_sub_epi64(vr_mantissa, vres_mantissa);	
	vr_lt_vres_exponent = _mm256_sub_epi64(vr_exponent, vres_exponent);
	
	k = _mm256_and_si256(vr_lt_vres_mantissa, vr_lt_vres_exponent);
	k = _mm256_or_si256(k, vres_eq_1);

	_mm256_store_si256((__m256i *)(b), k);
}

/* make sure that Pr(rerun the PRG)<=2^(-256) */
static void uniform_sampler_b1_k(const unsigned char *r, __m256i *y1, __m256i *y2)
{
	uint64_t sample[8] __attribute__ ((aligned (32)));
	uint32_t i = 0, j = 0;
	uint64_t x;
	
	while (j < 8)
	{
		do
		{	/* we ignore the low probability of rerunning the PRG 
			 * change the loading size i.e. uint16_t according to your UNIFORM_SIZE */
			x = load_16(r + UNIFORM_SIZE_B1_K * (i++));
		} while (1 ^ ((x - UNIFORM_Q_B1_K) >> 63));

		x = x - ((((x * BARRETT_FACTOR_B1_K) >> BARRETT_BITSHIFT_B1_K) + 1) * B1_K);
		x = x + (x >> 63) * B1_K;
		
		sample[j++] = x;
	}
	
	*y1 = _mm256_load_si256((__m256i *)(sample));
	*y2 = _mm256_load_si256((__m256i *)(sample + 4));
}

static void uniform_sampler_b2_k(const unsigned char *r, __m256i *y1, __m256i *y2)
{
	uint64_t sample[8] __attribute__ ((aligned (32)));
	uint32_t i = 0, j = 0;
	uint64_t x;
	
	while (j < 8)
	{
		do
		{	/* we ignore the low probability of rerunning the PRG 
			 * change the loading size i.e. uint16_t according to your UNIFORM_SIZE */
			x = load_16(r + UNIFORM_SIZE_B2_K * (i++));
		} while (1 ^ ((x - UNIFORM_Q_B2_K) >> 63));

		x = x - ((((x * BARRETT_FACTOR_B2_K) >> BARRETT_BITSHIFT_B2_K) + 1) * B2_K);
		x = x + (x >> 63) * B2_K;
		
		sample[j++] = x;
	}
	
	*y1 = _mm256_load_si256((__m256i *)(sample));
	*y2 = _mm256_load_si256((__m256i *)(sample + 4));
}

static void uniform_sampler_b3_k(const unsigned char *r, __m256i *y1, __m256i *y2)
{
	uint64_t sample[8] __attribute__ ((aligned (32)));
	uint32_t i = 0, j = 0;
	uint64_t x;
	
	while (j < 8)
	{
		do
		{	/* we ignore the low probability of rerunning the PRG 
			 * change the loading size i.e. uint16_t according to your UNIFORM_SIZE */
			x = load_24(r + UNIFORM_SIZE_B3_K * (i++));
		} while (1 ^ ((x - UNIFORM_Q_B3_K) >> 63));

		x = x - ((((x * BARRETT_FACTOR_B3_K) >> BARRETT_BITSHIFT_B3_K) + 1) * B3_K);
		x = x + (x >> 63) * B3_K;
		
		sample[j++] = x;
	}
	
	*y1 = _mm256_load_si256((__m256i *)(sample));
	*y2 = _mm256_load_si256((__m256i *)(sample + 4));
}

/* binary sampling algorithm 
 * we compute 8 samples every time by using the AVX2, 
 * then do the rejection */
static inline void gaussian_sampler(int64_t *sample, const uint32_t slen, const __m256i v_k, const __m256d v_k_2_inv, const uint64_t uniform_r_size, void (*uniform_sampler)(const unsigned char *, __m256i *, __m256i *))
{
	__m256i v_x, v_y1, v_y2, v_z, v_b_in;
	uint64_t z[8] __attribute__ ((aligned (32)));
	uint64_t b[8] __attribute__ ((aligned (32)));
	
	unsigned char r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_R_SIZE_MAX + 1] __attribute__ ((aligned (32)));
	unsigned char *r1;
	
	uint32_t i = 8, j = 0;
	uint64_t k;
	
	while (j < slen)
	{
		do
		{
			if (i == 8)
			{
				/* x<--D_sigma_0, y<--U([0,k-1]), z=kx+y */
				fastrandombytes_prv(r, 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + uniform_r_size + 1);
				
				uniform_sampler(r + 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE), &v_y1, &v_y2);
				
				r1 = r;
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, v_k);
				v_z = _mm256_add_epi64(v_x, v_y1);
				_mm256_store_si256((__m256i *)(z), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y1);
				bernoulli_sampler(b, v_b_in, r1 + BASE_TABLE_SIZE, v_k_2_inv);
				
				r1 = r + (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE);
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, v_k);
				v_z = _mm256_add_epi64(v_x, v_y2);
				_mm256_store_si256((__m256i *)(z + 4), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y2);
				bernoulli_sampler(b + 4, v_b_in, r1 + BASE_TABLE_SIZE, v_k_2_inv);

				i = 0;
			}
			
			k = (r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + uniform_r_size] >> i) & 0x1;
			i++;			
		} while (1 ^ ((b[i - 1] & ((z[i - 1] | -z[i - 1]) | (k | -k))) >> 63)); /* rejection condition: b=0 or ((b=1) && (z=0) && (k=0)) */
		
		sample[j++] = z[i - 1] * (1 ^ ((-k) & 0xfffffffffffffffe)); /* sample=z*(-1)^k */
	}
}

static inline void gaussian_b(POLY_R *out, const uint64_t m, const __m256i v_k, const __m256d v_k_2_inv, const uint64_t uniform_r_size, void (*uniform_sampler)(const unsigned char *, __m256i *, __m256i *))
{
	uint64_t i;
	
	for (i = 0; i < m; i++)
	{
		gaussian_sampler((int64_t *)((out + i)->poly), D, v_k, v_k_2_inv, uniform_r_size, uniform_sampler);
	}
}

void sample_b1(POLY_R *out)
{
	gaussian_b(out, 2 * L * NHAT * K, V_B1_K, V_B1_K_2_INV, UNIFORM_REJ_B1_K * UNIFORM_SIZE_B1_K, uniform_sampler_b1_k);
}

void sample_b2(POLY_R *out)
{
	gaussian_b(out, L, V_B2_K, V_B2_K_2_INV, UNIFORM_REJ_B2_K * UNIFORM_SIZE_B2_K, uniform_sampler_b2_k);
}

void sample_b3(POLY_R *out)
{
	gaussian_b(out, M, V_B3_K, V_B3_K_2_INV, UNIFORM_REJ_B3_K * UNIFORM_SIZE_B3_K, uniform_sampler_b3_k);
}

/* Bernoulli sampler with bias exp((-2 * <z, c> + ||c||^2) / (2 * \sigma^2)) / mu(\phi) in Rej/RejOp */
static inline uint64_t rej_bernoulli_sampler(uint64_t x, const unsigned char *r, const __m128d v_k_2_inv, const __m128d v_mu)
{
	/* coefficients of the exp evaluation polynomial */
	static const uint64_t EXP_COFF[] = {0x3e833b70ffa2c5d4,
										0x3eb4a480fda7e6e1,
										0x3ef01b254493363f,
										0x3f242e0e0aa273cc,
										0x3f55d8a2334ed31b,
										0x3f83b2aa56db0f1a,
										0x3fac6b08e11fc57e,
										0x3fcebfbdff556072,
										0x3fe62e42fefa7fe6,
										0x3ff0000000000000};

	__m128d vx, vx_1, vx_2, vsum;
	__m128i vt, vres;

	uint64_t res, res_mantissa, res_exponent;
	uint64_t r_mantissa, r_exponent;
	uint64_t r1, r2;
	
	/* 2^x=2^(floor(x)+a)=2^(floor(x))*2^a, where a is in [0,1]
	 * we only evaluate 2^a by using a polynomial */	
	vx = _mm_cvtsi64_sd(_mm_setzero_pd(), x);
	vx = _mm_mul_sd(vx, v_k_2_inv);
	vx = _mm_sub_sd(vx, v_mu);
	
	vx_1 = _mm_floor_pd(vx);
	vt = _mm_cvtpd_epi32(vx_1);
	vt = _mm_slli_epi64(vt, 52);
	
	/* evaluate 2^a */
	vx_2 = _mm_sub_sd(vx, vx_1);
	vsum = _mm_add_sd(_mm_mul_sd(_mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[0])), vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[1])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[2])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[3])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[4])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[5])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[6])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[7])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[8])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx_2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[9])));
	
	/* combine to compute 2^x */
	vres = _mm_add_epi64(vt, _mm_castpd_si128(vsum));
	_mm_storel_epi64((__m128i *)(&res), vres);
	
	res_mantissa = (res & EXP_MANTISSA_MASK) | (1LL << EXP_MANTISSA_PRECISION);
	res_exponent = REJ_R_EXPONENT_L - 1023 + 1 + (res >> EXP_MANTISSA_PRECISION); 
	
	r1 = *((uint64_t *)r);
	r2 = (uint64_t)(r[8]);
	
	r_mantissa = r1 & R_MANTISSA_MASK;
	r_exponent = (r1 >> R_MANTISSA_PRECISION) | (r2 << (64 - R_MANTISSA_PRECISION));
	
	/* ~((res == 1.0) || ((r_mantissa < res_mantissa) && (r_exponent < (1 << res_exponent)))) */
	return (((res - DOUBLE_ONE) | (DOUBLE_ONE - res)) & ((1LL << 63) ^ ((r_mantissa - res_mantissa) & (r_exponent - (1LL << res_exponent))))) >> 63;
}

static inline uint64_t rej(const POLY_R *z, const POLY_R *c, const uint64_t slen, const __m128d v_k_2_inv, const __m128d v_mu, const uint64_t sigma_6)
{
	unsigned char r[REJ_BERNOULLI_ENTRY_SIZE];
	
	uint64_t ret = 0, x = 0;
	
	uint64_t i, j;
	
	for (i = 0; i < slen; i++)
	{
		for (j = 0; j < D; j++)
		{
			/* check z norm */
			ret |= ct_lt(sigma_6, z[i].poly[j]) | ct_lt(z[i].poly[j], -sigma_6);
			
			/* -2 * <z, c> + ||c||^2 */
			x += (int64_t)(c[i].poly[j]) * (int64_t)(c[i].poly[j]) - 2 * (int64_t)(z[i].poly[j]) * (int64_t)(c[i].poly[j]);
		}
	}
	
	fastrandombytes_prv(r, REJ_BERNOULLI_ENTRY_SIZE);
	ret |= rej_bernoulli_sampler(x, r, v_k_2_inv, v_mu);
	
	return ret;
}

uint64_t rej_op(const POLY_R *z, const POLY_R *z_a, const POLY_R *z_J, const POLY_R *x_r, const POLY_R *x_b, const POLY_R *x_re)
{
	unsigned char r[REJ_BERNOULLI_ENTRY_SIZE];
	
	uint64_t ret = 0, x = 0;
	
	uint64_t i, j;
	
	/* <z, c> ?< 0 */
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < D; j++)
		{
			x += (int64_t)(z[i].poly[j]) * (int64_t)(x_r[i].poly[j]);
		}
	}
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < D; j++)
		{
			x += (int64_t)(z_a[i].poly[j]) * (int64_t)(x_b[i].poly[j]);
		}
	}
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < D; j++)
		{
			x += (int64_t)(z_J[i].poly[j]) * (int64_t)(x_re[i].poly[j]);
		}
	}
	
	if (x >> 63)
	{
		return 1;
	}
	
	x *= -2;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < D; j++)
		{
			ret |= ct_lt(B3_6, z[i].poly[j]) | ct_lt(z[i].poly[j], -B3_6);

			x += (int64_t)(x_r[i].poly[j]) * (int64_t)(x_r[i].poly[j]);
		}
	}
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < D; j++)
		{
			ret |= ct_lt(B3_6, z_a[i].poly[j]) | ct_lt(z_a[i].poly[j], -B3_6);

			x += (int64_t)(x_b[i].poly[j]) * (int64_t)(x_b[i].poly[j]);
		}
	}
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < D; j++)
		{
			ret |= ct_lt(B3_6, z_J[i].poly[j]) | ct_lt(z_J[i].poly[j], -B3_6);

			x += (int64_t)(x_re[i].poly[j]) * (int64_t)(x_re[i].poly[j]);
		}
	}
	
	fastrandombytes_prv(r, REJ_BERNOULLI_ENTRY_SIZE);
	ret |= rej_bernoulli_sampler(x, r, REJ_V_B3_K_2_INV, V_MU_B3);
	
	return ret;
}

uint64_t rej_f(const POLY_R *z, const POLY_R *c)
{
	return rej(z, c, 2 * L * NHAT * K, REJ_V_B1_K_2_INV, V_MU_B1, B1_6);
}

uint64_t rej_fJ(const POLY_R *z, const POLY_R *c)
{
	return rej(z, c, L, REJ_V_B2_K_2_INV, V_MU_B2, B2_6);
}
