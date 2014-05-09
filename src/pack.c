/*
 * CPQREF/pack.c
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

#include "params.h"
#include "pack.h"
#include "pqerror.h"


int
get_public_key_ptrs(
    PQ_PARAM_SET        **params,
    int64_t             **h,
    unsigned char       **digest,
    const size_t        public_key_len,
    const unsigned char *public_key_blob)
{
  PQ_PARAM_SET *P;
  int64_t *lh;
  unsigned char *ldigest;
  size_t expected_size;

  if(!public_key_blob || public_key_len < OID_BYTES)
  {
    return PQNTRU_ERROR;
  }

  if(!(P = pq_get_param_set_by_oid((const uint8_t*)public_key_blob)))
  {
    return PQNTRU_ERROR;
  }

  expected_size = PUB_KEY_PACKED_BYTES(P);
  if(public_key_len != expected_size)
  {
    return PQNTRU_ERROR;
  }

  lh = (int64_t*)(public_key_blob + OID_BYTES);
  ldigest = (unsigned char*)(lh + P->N);

  if(params)
  {
    *params = P;
  }

  if(h)
  {
    *h = lh;
  }

  if(digest)
  {
    *digest = ldigest;
  }

  return PQNTRU_OK;
}


int
get_private_key_ptrs(
    PQ_PARAM_SET        **params,
    uint16_t            **fi,
    uint16_t            **gi,
    int8_t              **ginv,
    const size_t        private_key_len,
    const unsigned char *private_key_blob)
{
  size_t df;
  PQ_PARAM_SET *P;
  uint16_t *lfi;
  uint16_t *lgi;
  int8_t *lginv;
  size_t expected_size;

  if(!private_key_blob || private_key_len < OID_BYTES)
  {
    return PQNTRU_ERROR;
  }

  if(!(P = pq_get_param_set_by_oid((const uint8_t*)private_key_blob)))
  {
    return PQNTRU_ERROR;
  }

  expected_size = PRIV_KEY_PACKED_BYTES(P);
  if(private_key_len != expected_size)
  {
    return PQNTRU_ERROR;
  }

  df = 2*(P->d1 + P->d2 + P->d3);

  lfi = (uint16_t*)(private_key_blob + OID_BYTES);
  lgi = lfi + df;
  lginv = (int8_t*)(lgi + df);

  if(params)
  {
    *params = P;
  }

  if(fi)
  {
    *fi = lfi;
  }

  if(gi)
  {
    *gi = lgi;
  }

  if(ginv)
  {
    *ginv = lginv;
  }

  return PQNTRU_OK;
}

