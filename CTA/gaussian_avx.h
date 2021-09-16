#ifndef _GAUSSIAN_AVX_H
#define _GAUSSIAN_AVX_H

#include <stdint.h>
#include "poly_q.h"
#include "param.h"

void sample_b1(POLY_R *out);
void sample_b2(POLY_R *out);
void sample_b3(POLY_R *out);

uint64_t rej_f(const POLY_R *z, const POLY_R *c);
uint64_t rej_fJ(const POLY_R *z, const POLY_R *c);

//uint64_t rej(const POLY_R *z, const POLY_R *c, const uint64_t slen, const __m128d v_k_2_inv, const __m128d v_mu, const uint64_t sigma_6);
uint64_t rej_op(const POLY_R *z, const POLY_R *z_a, const POLY_R *z_J, const POLY_R *x_r, const POLY_R *x_b, const POLY_R *x_re);

#endif
