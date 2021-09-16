#ifndef _POLY_RED_H
#define _POLY_RED_H

#include <stdint.h>
#include "poly_param.h"

#define RED_SHORT_SHIFT_1_Q 32LL
#define RED_SHORT_CONST 0LL
#define RED_SHORT_MASK_Q ((1LL << RED_SHORT_SHIFT_1_Q) - 1)

#define RED_SHORT_SHIFT_1_QHAT 6LL
#define RED_SHORT_MASK_QHAT ((1LL << RED_SHORT_SHIFT_1_QHAT) - 1)

static inline uint64_t con_sub(const uint64_t x, const uint64_t q)
{
	return x - ((-(1 ^ ((x - q) >> 63))) & q);
}

static inline uint64_t con_add(const uint64_t x, const uint64_t q)
{
	return x + ((-(x >> 63)) & q);
}

static inline uint64_t red_short_q(uint64_t t)
{
	return con_add((t & RED_SHORT_MASK_Q) - ((t >> RED_SHORT_SHIFT_1_Q) * RED_SHORT_CONST), Q);
}

static inline uint64_t red_short_qhat(uint64_t t)
{
	uint64_t x, y;
	
	x = t >> RED_SHORT_SHIFT_1_QHAT;
	y = t & RED_SHORT_MASK_QHAT;
	
	return con_sub(y - x, QHAT);
}

#define MONTGOMERY_FACTOR_Q 3049919425LL
#define MONTGOMERY_SHIFT_Q 32LL
#define MONTGOMERY_MASK_Q ((1LL << MONTGOMERY_SHIFT_Q) - 1)
#define MONTGOMERY_CONVERT_FACTOR_Q 919681LL
#define MONTGOMERY_32_Q 959LL

//not check
#define MONTGOMERY_FACTOR_QHAT 337657985LL
#define MONTGOMERY_SHIFT_QHAT 30LL
#define MONTGOMERY_MASK_QHAT ((1LL << MONTGOMERY_SHIFT_QHAT) - 1)
#define MONTGOMERY_CONVERT_FACTOR_QHAT 8824684028LL
#define MONTGOMERY_SEP_QHAT_SHIFT 31LL
#define MONTGOMERY_SEP_QHAT_MASK ((1LL << MONTGOMERY_SEP_QHAT_SHIFT) - 1)
#define MONTGOMERY_SEP_QHAT_HI 4294967296LL

/* Montgomery reduction
 * Input: x < Q*R, where R=2^k
 * Output: m = x*R^{-1} % Q
 * 
 * b = Q^{-1} % R
 * t = ((x % R)*b) % R
 * m = x / R - t * Q / R */
static inline uint64_t montgomery(uint64_t a, uint64_t b, const uint64_t q, const uint64_t montgomery_factor, const uint64_t montgomery_shift, const uint64_t montgomery_mask)
{
	uint64_t t, x, y;
	
	t = a * b;
	x = t & montgomery_mask;
	y = (x * montgomery_factor) & montgomery_mask;

	return con_add((t >> montgomery_shift) - ((y * q) >> montgomery_shift), q);
}
 
static inline uint64_t montgomery_q(uint64_t a, uint64_t b)
{
	return montgomery(a, b, Q, MONTGOMERY_FACTOR_Q, MONTGOMERY_SHIFT_Q, MONTGOMERY_MASK_Q);
}

static inline uint64_t montgomery_qhat(uint64_t a, uint64_t b)
{
	uint64_t a0, a1, b0, b1;
	uint64_t x0, x1, x2;
	
	a0 = a & MONTGOMERY_SEP_QHAT_MASK;
	a1 = a >> MONTGOMERY_SEP_QHAT_SHIFT;
	b0 = b & MONTGOMERY_SEP_QHAT_MASK;
	b1 = b >> MONTGOMERY_SEP_QHAT_SHIFT;
	
	x0 = montgomery(a0, b0, QHAT, MONTGOMERY_FACTOR_QHAT, MONTGOMERY_SHIFT_QHAT, MONTGOMERY_MASK_QHAT);
	x1 = a0 * b1 + a1 * b0;
	x2 = (a1 * b1 + (x1 >> MONTGOMERY_SEP_QHAT_SHIFT)) * MONTGOMERY_SEP_QHAT_HI;
	x1 = montgomery(x1 & MONTGOMERY_SEP_QHAT_MASK, 1LL << MONTGOMERY_SEP_QHAT_SHIFT, QHAT, MONTGOMERY_FACTOR_QHAT, MONTGOMERY_SHIFT_QHAT, MONTGOMERY_MASK_QHAT);
	
	return red_short_qhat(x0 + x1 + x2);
}

/* Barrett reduction
 * Input: x < 2^k
 * Output m = x % Q in [0, 2Q)
 * 
 * b = floor(2^k/Q)
 * t = floor((x * b) / 2^k), where t is an estimation of x / Q
 * m = x - t * Q */
static inline uint64_t barrett(const uint64_t x, const uint64_t q, const uint64_t barrett_factor, const uint64_t barrett_shift)
{
	return con_sub(x - ((x * barrett_factor) >> barrett_shift) * q, q);
}

#endif

