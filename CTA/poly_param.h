#ifndef _POLY_PARAM_H
#define _POLY_PARAM_H

#define R 32LL

#define D 128LL /* ring dimension */

/* modulus */
#define Q 4294966337LL
#define Q_BYTE 4LL
#define STORE_Q store_32

#define QHAT 61LL
#define QHAT_BYTE 1LL
#define STORE_QHAT store_8

#define SUBRING_SIZE (D / R)

#endif
