/*
 * CPQREF/pqntrusign.c
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


#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "crypto_hash_sha512.h"

#include "shred.h"
#include "params.h"
#include "pack.h"
#include "pol.h"

#include "pqerror.h"

#define STATS 1

int g_loop = 0;
int g_max = 0;
int g_max2 = 0;

static int
challenge(
    int8_t              *sp,
    int8_t              *tp,
    const size_t        public_key_len,
    const unsigned char *public_key_blob,
    const size_t        msg_len,
    const unsigned char *msg)
{
  uint16_t      i;
  uint16_t      j;
  uint16_t      m;

  uint8_t       r;
  int8_t        *poly;

  PQ_PARAM_SET  *P;

  unsigned char *digest;

  unsigned char input[2*HASH_BYTES];
  unsigned char pool[HASH_BYTES];

  int result;

  result = get_public_key_ptrs(&P, NULL, &digest,
                               public_key_len, public_key_blob);
  if(PQNTRU_ERROR == result)
  {
    return PQNTRU_ERROR;
  }


  /* pool = hash(hash(msg) || hash(public key)) */
  crypto_hash_sha512(input, msg, msg_len);
  memcpy(input+HASH_BYTES, digest, HASH_BYTES);
  crypto_hash_sha512(pool, input, 2*HASH_BYTES);

  j = 0;
  poly = sp;
  for(m = 0; m < 2; m++)
  {
    i = 0;
    r = 0;
    while(i < P->N)
    {
      if(HASH_BYTES == j)
      {
        memcpy(input, pool, HASH_BYTES);
        crypto_hash_sha512(pool, input, HASH_BYTES);
        j = 0;
      }
      if(0 == r)
      {
        r = (uint8_t) pool[j++];
      }
      switch(r & 3)
      {
        case 1: poly[i] = -1; break;
        case 2: poly[i] =  0; break;
        case 3: poly[i] =  1; break;
        default:  r >>= 2; continue;
      }
      r >>= 2;
      i++;
    }

    poly = tp;
  }

  return 0;
}



int
pq_sign(
    int64_t             *sig,
    const size_t        private_key_len,
    const unsigned char *private_key_blob,
    const size_t        public_key_len,
    const unsigned char *public_key_blob,
    const size_t        msg_len,
    const unsigned char *msg)
{

  uint16_t i;
  int error = 0;

  uint16_t  N;
  uint16_t  padN;
  int64_t   q;
  int8_t    p;
  uint16_t  d1;
  uint16_t  d2;
  uint16_t  d3;
  int64_t   m;

  size_t        scratch_len;
  unsigned char *scratch;

  int8_t  *sp;
  int8_t  *tp;
  int64_t *s0;
  int64_t *t0;
  int64_t *ginv64;
  int64_t *a;
  int64_t *tmpx2;

  PQ_PARAM_SET  *P;
  uint16_t      *fi;
  uint16_t      *gi;
  int8_t        *ginv;
  int64_t       *h;
  unsigned char *digest;

  int result = PQNTRU_OK;

  result = get_private_key_ptrs(&P, &fi, &gi, &ginv,
                                private_key_len, private_key_blob);
  if(result == PQNTRU_ERROR)
  {
    return PQNTRU_ERROR;
  }

  result = get_public_key_ptrs(NULL, &h, &digest,
                               public_key_len, public_key_blob);
  if(result == PQNTRU_ERROR)
  {
    return PQNTRU_ERROR;
  }

  N = P->N;
  padN = P->padded_N;
  q = P->q;
  p = P->p;
  d1 = P->d1;
  d2 = P->d2;
  d3 = P->d3;

  scratch_len = 6 * POLYNOMIAL_BYTES(P) + 2 * N;
  if(!(scratch = malloc(scratch_len)))
  {
    return PQNTRU_ERROR;
  }
  memset(scratch, 0, scratch_len);

  s0     = (int64_t*)(scratch);
  t0     = (int64_t*)(scratch + POLYNOMIAL_BYTES(P));
  ginv64 = (int64_t*)(scratch + 2*POLYNOMIAL_BYTES(P));
  /* a is treated as 3 polynomials, aliases tmpx2 */
  a      = (int64_t*)(scratch + 3*POLYNOMIAL_BYTES(P));
  tmpx2  = (int64_t*)(scratch + 4*POLYNOMIAL_BYTES(P));
  sp     =  (int8_t*)(scratch + 6*POLYNOMIAL_BYTES(P));
  tp     =  (int8_t*)(scratch + 6*POLYNOMIAL_BYTES(P) + N);

  for(i=0; i<N; i++)
  {
    ginv64[i] = (int64_t) ginv[i];
  }

  challenge(sp, tp,
            public_key_len, public_key_blob,
            msg_len, msg);

#if STATS
  int loop = 0;
#endif
  do
  {
    error = 0;

    /* Choose random s0 satisfying s0 = sp (mod p) */
    pol_unidrnd_pZ(s0, N, q, p);
    for(i=0; i<N; i++)
    {
      s0[i] += sp[i];
    }

    /* Load h into a zero padded polynomial */
    memcpy(t0, h, N*sizeof(int64_t));

    /* t0 = h*s0 */
    pol_mul_coefficients(t0, t0, s0, N, padN, q, a);

    /* t0 = tp - (s0*h) */
    for(i=0; i<N; i++)
    {
      t0[i] *= -1;
      t0[i] += tp[i];
    }

    /* a = ginv * (tp - t0) (mod p) */
    pol_mul_coefficients(a, t0, ginv64, N, padN, p, a);

    /* tmpx2 = a * F = (a * (f-1)/p) */
    pol_mul_product(tmpx2, a, d1, d2, d3, fi, N, tmpx2);
    for(i=0; i<N; i++)
    {
      m = p * (a[i] + tmpx2[i]);
      error |= (m > P->B_s) || (-m > P->B_s);

      /* s0 = s0 + p*(a + tmpx2) = s0 + a*f */
      s0[i] += m;

      error |= (cmod(s0[i], p) - sp[i]); /* Not necessary to check this */
      error |= (s0[i] > P->norm_bound_s) || (-s0[i] > P->norm_bound_s);
#if STATS
      if(m > g_max)
      {
        g_max = m;
      }
      else if(-m > g_max)
      {
        g_max =  -m;
      }
#endif
    }

    /* tmpx2 = a * G = (a * (g - 1)) */
    pol_mul_product(tmpx2, a, d1, d2, d3, gi, N, tmpx2);
    for(i=0; i<N; i++)
    {
      m = (a[i] + tmpx2[i]);
      error |= (m > P->B_t) || (-m > P->B_t);

      /* t0 = (a + tmpx2) - t0 + tp = a*g - tp + s0*h + tp = s0*h + a*g */
      t0[i] = m - t0[i] + tp[i];
      error |= (cmod(t0[i], p) - tp[i]); /* Not necessary to check this */
      error |= (t0[i] > P->norm_bound_t) || (-t0[i] > P->norm_bound_t) ;
#if STATS
      if(m > g_max2)
      {
        g_max2 = m;
      }
      else if(-m > g_max2)
      {
        g_max2 = -m;
      }
#endif
    }
#if STATS
    loop++;
#endif
  } while(0 != error);

#if STATS
  g_loop += loop;
#endif
  memcpy(sig, s0, N * sizeof(int64_t));

  shred(scratch, scratch_len);
  free(scratch);

  return PQNTRU_OK;
}


int
pq_verify(
    int64_t             *sig,
    const size_t        public_key_len,
    const unsigned char *public_key_blob,
    const size_t        msg_len,
    const unsigned char *msg)
{
  uint16_t      i;

  PQ_PARAM_SET  *P;
  size_t        scratch_len;
  unsigned char *scratch;
  int8_t        *sp;
  int8_t        *tp;
  int64_t       *h;
  int64_t       *lh;
  int64_t       *lsig;
  int64_t       *tmpx3;

  uint16_t      N;
  uint16_t      padN;
  int64_t       q;
  int8_t        p;
  unsigned char *digest;
  uint16_t      error = 0;
  int           result;


  result = get_public_key_ptrs(&P, &h, &digest,
                               public_key_len, public_key_blob);

  if(PQNTRU_ERROR == result)
  {
    return PQNTRU_ERROR;
  }

  N = P->N;
  padN = P->padded_N;
  p = P->p;
  q = P->q;

  scratch_len = 2 * N + 5 * POLYNOMIAL_BYTES(P);
  if(!(scratch = malloc(scratch_len)))
  {
    return PQNTRU_ERROR;
  }
  memset(scratch, 0, scratch_len);

  lsig  = (int64_t*)scratch;
  lh    = (int64_t*)(scratch + POLYNOMIAL_BYTES(P));
  tmpx3 = (int64_t*)(scratch + 2*POLYNOMIAL_BYTES(P));
  sp    =  (int8_t*)(scratch + 5*POLYNOMIAL_BYTES(P));
  tp    =  (int8_t*)(scratch + 5*POLYNOMIAL_BYTES(P) + N);

  memcpy(lh, h, N*sizeof(int64_t));
  memcpy(lsig, sig, N*sizeof(int64_t));

  challenge(sp, tp,
            public_key_len, public_key_blob,
            msg_len, msg);

  pol_mul_coefficients(tmpx3, lsig, lh, N, padN, q, tmpx3);

  for(i=0; i<N; i++)
  {
    error |= (cmod(sig[i], p) - sp[i]);
    error |= (sig[i] > P->norm_bound_s) || (-sig[i] > P->norm_bound_s);

    error |= (cmod(tmpx3[i], p) - tp[i]);
    error |= (tmpx3[i] > P->norm_bound_t) || (-tmpx3[i] > P->norm_bound_t) ;
  }

  free(scratch);

  if(0 == error)
  {
    return PQNTRU_OK;
  }
  return PQNTRU_ERROR;
}


int
pq_gen_key(
    PQ_PARAM_SET  *P,
    size_t        *privkey_blob_len,
    unsigned char *privkey_blob,
    size_t        *pubkey_blob_len,
    unsigned char *pubkey_blob)
{
  uint16_t      i;
  uint16_t      m;

  uint16_t      N;
  uint16_t      padN;
  int64_t       q;
  int8_t        p;
  uint16_t      d1;
  uint16_t      d2;
  uint16_t      d3;

  size_t        private_key_blob_len;
  size_t        public_key_blob_len;

  uint16_t      *fi;
  uint16_t      *gi;
  int8_t        *ginv;
  int64_t       *h;
  unsigned char *digest;

  size_t        scratch_len;
  unsigned char *scratch;
  int64_t       *a1;
  int64_t       *a2;
  int64_t       *tmpx3;

  if(!P || !privkey_blob_len || !pubkey_blob_len)
  {
    return PQNTRU_ERROR;
  }

  N = P->N;
  padN = P->padded_N;
  q = P->q;
  p = P->p;
  d1 = P->d1;
  d2 = P->d2;
  d3 = P->d3;

  /* TODO: Standardize packed key formats */

  private_key_blob_len = PRIV_KEY_PACKED_BYTES(P);
  public_key_blob_len = PUB_KEY_PACKED_BYTES(P);

  if(!privkey_blob || !pubkey_blob)
  {
    if(!privkey_blob && privkey_blob_len != NULL)
    {
      *privkey_blob_len = private_key_blob_len;
    }
    if(!pubkey_blob && pubkey_blob_len != NULL)
    {
      *pubkey_blob_len = public_key_blob_len;
    }
    return PQNTRU_OK;
  }

  if((*privkey_blob_len != private_key_blob_len)
      || (*pubkey_blob_len != public_key_blob_len))
  {
    return PQNTRU_ERROR;
  }

  memcpy(privkey_blob, P->OID, OID_BYTES);
  memcpy(pubkey_blob, P->OID, OID_BYTES);

  get_private_key_ptrs(NULL, &fi, &gi, &ginv,
                       private_key_blob_len, privkey_blob);
  get_public_key_ptrs(NULL, &h, &digest,
                      public_key_blob_len, pubkey_blob);

  scratch_len = 5 * POLYNOMIAL_BYTES(P);

  if(!(scratch = malloc(scratch_len)))
  {
    return PQNTRU_ERROR;
  }
  a1 = (int64_t*)scratch;
  a2 = (int64_t*)(scratch + POLYNOMIAL_BYTES(P));
  tmpx3 = (int64_t*)(scratch + 2*POLYNOMIAL_BYTES(P));;
  memset(scratch, 0, scratch_len);


  /* Find invertible pf mod q */
  /* TODO: Better sampling of product form keys
   *       Try to avoid keys with f(1) = 0
   */
  do
  {
    pol_gen_product(fi, d1, d2, d3, N);

    /* f = p * (1 + product form poly) */
    memset(a1, 0, POLYNOMIAL_BYTES(P));
    a1[0] = p;

    pol_mul_product(a1, a1, d1, d2, d3, fi, N, tmpx3);
    a1[0] += p;

  } while(PQNTRU_ERROR == pol_inv_mod2(a2, a1, N));

  /* Lift from (Z/2Z)[X]/(X^N - 1) to (Z/qZ)[X]/(X^N -1) */
  for (m = 0; m < 5; ++m)   /* assumes 2^16 < q <= 2^32 */
  {
    /* a^-1 = a^-1 * (2 - a * a^-1) mod q */

    pol_mul_product(a1, a2, d1, d2, d3, fi, N, tmpx3);

    for (i = 0; i < N; ++i)
    {
      a1[i] = -p*(a1[i] + a2[i]);
    }

    a1[0] = a1[0] + 2;
    pol_mul_coefficients(a2, a2, a1, N, padN, q, tmpx3);
  }


  /* Find invertible g mod p */
  do
  {
    /* Generate product form g,
     * then expand it to find inverse mod p
     */
    pol_gen_product(gi, d1, d2, d3, N);

    memset(a1, 0, POLYNOMIAL_BYTES(P));
    a1[0] = 1;

    pol_mul_product(a1, a1, d1, d2, d3, gi, N, tmpx3);
    a1[0] += 1;

  } while(PQNTRU_ERROR == pol_inv_modp(tmpx3, a1, N, p));

  for(i=0; i<N; i++)
  {
    ginv[i] = (int8_t)(tmpx3[i] & 0xff);
  }

  /* Calculate h = g/f mod q */
  pol_mul_product(h, a2, d1, d2, d3, gi, N, tmpx3);
  for(i=0; i<N; i++)
  {
    h[i] = cmod(h[i] + a2[i], q);
  }

  crypto_hash_sha512(digest, pubkey_blob, N*sizeof(int64_t));

  shred(scratch, scratch_len);
  free(scratch);

  return PQNTRU_OK;
}


