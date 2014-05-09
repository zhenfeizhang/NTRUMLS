/*
 * CPQREF/pack.h
 *
 *  Copyright 2014 John M. Schanck
 *
 *  This file is part of CPQREF.
 *
 *  CPQREF is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CPQREF is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CPQREF.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CPQREF_PACK_H_
#define CPQREF_PACK_H_


#define HASH_BYTES 64

#define OID_BYTES ((sizeof(((PQ_PARAM_SET*)0)->OID)))

#define POLYNOMIAL_BYTES(P) ((P)->padded_N * sizeof(int64_t))

/* Private key := OID || f indices || g indices || g inverse */
#define PRIV_KEY_PACKED_BYTES(P) \
         (OID_BYTES \
         + 2*((P)->d1 + (P)->d2 + (P)->d3) * sizeof(uint16_t) \
         + 2*((P)->d1 + (P)->d2 + (P)->d3) * sizeof(uint16_t) \
         + (P)->N)

/* Public key := OID || h || hash(h) */
#define PUB_KEY_PACKED_BYTES(P)  \
        (OID_BYTES + (P)->N * sizeof(int64_t) + HASH_BYTES)


int
get_public_key_ptrs(
    PQ_PARAM_SET        **params,
    int64_t             **h,
    unsigned char       **digest,
    const size_t        public_key_len,
    const unsigned char *public_key_blob);

int
get_private_key_ptrs(
    PQ_PARAM_SET        **params,
    uint16_t            **fi,
    uint16_t            **gi,
    int8_t              **ginv,
    const size_t        public_key_len,
    const unsigned char *public_key_blob);

#endif
