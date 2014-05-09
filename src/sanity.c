/*
 * CPQREF/sanity.c
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
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "fastrandombytes.h"
#include "params.h"
#include "pqerror.h"
#include "pqntrusign.h"

#include "pack.h"
#include "pol.h"

#define TIMES 1000

static int
testKeyGen(PQ_PARAM_SET_ID id)
{
  uint16_t i;
  uint16_t j;

  PQ_PARAM_SET *P;
  size_t privkey_blob_len;
  size_t pubkey_blob_len;

  unsigned char *privkey_blob;
  unsigned char *pubkey_blob;

  unsigned char *scratch;
  size_t         scratch_len;

  uint16_t      *fi;
  uint16_t      *gi;
  int8_t        *ginv;
  int64_t       *h;

  if(!(P = pq_get_param_set_by_id(id)))
  {
    return -1;
  }

  scratch_len = 6*POLYNOMIAL_BYTES(P);
  scratch = malloc(scratch_len);
  int64_t *a1 = (int64_t*)scratch;
  int64_t *a2 = (int64_t*)(scratch + POLYNOMIAL_BYTES(P));
  int64_t *a3 = (int64_t*)(scratch + 2*POLYNOMIAL_BYTES(P));

  for(i=0; i<TIMES; i++)
  {
    memset(scratch, 0, 6 * POLYNOMIAL_BYTES(P));
    pq_gen_key(P, &privkey_blob_len, NULL, &pubkey_blob_len, NULL);

    privkey_blob = malloc(privkey_blob_len);
    pubkey_blob = malloc(pubkey_blob_len);

    if(PQNTRU_ERROR == pq_gen_key(P,
               &privkey_blob_len, privkey_blob,
               &pubkey_blob_len, pubkey_blob))
    {
      fprintf(stderr, "\t fail in keygen\n");
    }


    get_private_key_ptrs(NULL, &fi, &gi, &ginv, privkey_blob_len, privkey_blob);
    get_public_key_ptrs(NULL, &h, NULL, pubkey_blob_len, pubkey_blob);

    pol_mul_product(a1, h, P->d1, P->d2, P->d3, fi, P->N, a3);
    for(j=0; j<P->N; j++)
    {
      a1[j] = cmod(P->p * (h[j] + a1[j]), P->q);
    }

    for(j=0; j<P->N; j++)
    {
      a2[j] = (int64_t) ginv[j];
    }
    pol_mul_coefficients(a3, a1, a2, P->N, P->padded_N, P->p, a3);
    for(j=1; j<P->N; j++)
    {
      if(a3[j] != 0)
      {
        fprintf(stderr, "\t bad key");
        free(privkey_blob);
        free(pubkey_blob);
        return -1;
      }
    }

    free(privkey_blob);
    free(pubkey_blob);
  }

  return 0;
}


static int
testSet(PQ_PARAM_SET_ID id)
{
  uint16_t i;

  PQ_PARAM_SET *P;
  size_t privkey_blob_len;
  size_t pubkey_blob_len;
  unsigned char *privkey_blob;
  unsigned char *pubkey_blob;

  int64_t *sigs;
  uint16_t msg_len = 256;
  unsigned char *msg;

  int result = 0;

  if(!(P = pq_get_param_set_by_id(id)))
  {
    return -1;
  }
  fprintf(stderr, "Testing parameter set %s", P->name);
  fflush(stderr);

  pq_gen_key(P, &privkey_blob_len, NULL, &pubkey_blob_len, NULL);

  privkey_blob = malloc(TIMES * privkey_blob_len);
  pubkey_blob = malloc(TIMES * pubkey_blob_len);

  sigs = malloc(TIMES * P->N * sizeof(int64_t));
  msg = malloc(TIMES * msg_len * sizeof(int64_t));
  memset(msg, 0, TIMES*msg_len*sizeof(int64_t));


  for(i=0; i<TIMES; i++)
  {
    if(PQNTRU_ERROR == pq_gen_key(P,
               &privkey_blob_len, privkey_blob + (i*privkey_blob_len),
               &pubkey_blob_len, pubkey_blob + (i*pubkey_blob_len)))
    {
      result = -1;
      fprintf(stderr, "\t fail in keygen\n");
      goto exit;
    }
  }

  for(i=0; i<TIMES; i++)
  {
    fastrandombytes(msg+(i*msg_len), msg_len);
    if(PQNTRU_ERROR == pq_sign(sigs + (i * P->N),
                               privkey_blob_len, privkey_blob + (i*privkey_blob_len),
                               pubkey_blob_len, pubkey_blob + (i*pubkey_blob_len),
                               msg_len, msg + (i*msg_len)))
    {
      result = -1;
      fprintf(stderr, "\t fail in sign\n");
      goto exit;
    }
  }

  for(i=0; i<TIMES; i++)
  {
    if(PQNTRU_ERROR == pq_verify(sigs + (i * P->N),
                                 pubkey_blob_len, pubkey_blob + (i*pubkey_blob_len),
                                 msg_len, msg + (i*msg_len)))
    {
      result = -1;
      fprintf(stderr, "\t fail in verify\n");
      goto exit;
    }
  }

  fprintf(stderr, "\t good\n");

exit:
  free(msg);
  free(sigs);
  free(privkey_blob);
  free(pubkey_blob);
  return result;
}


int
main(int argc, char **argv)
{
  uint16_t i;
  PQ_PARAM_SET_ID plist[] = {DRAFT_401, DRAFT_439, DRAFT_593, DRAFT_743};
  size_t numParams = sizeof(plist)/sizeof(PQ_PARAM_SET_ID);

  for(i = 0; i<numParams; i++)
  {
    testKeyGen(plist[i]);
    testSet(plist[i]);
  }

  rng_cleanup();

  exit(EXIT_SUCCESS);
}
