/*
 * CPQREF/pol.c
 *
 *  Copyright 2013 John M. Schanck
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
#include <string.h>

#include "fastrandombytes.h"
#include "shred.h"
#include "pol.h"
#include "pqerror.h"



static int64_t
zz_gcd(
    const int64_t a,
    const int64_t b,
    int64_t       *u_ptr,
    int64_t       *v_ptr)
{
  int64_t d = a;
  int64_t u = 1;
  int64_t v = 0;
  int64_t v1, v3, t1, t3;
  if(b != 0) {
    v1 = 0;
    v3 = b;
    do {
      t1 = d / v3;
      t3 = d % v3;
      t1 = u - (t1*v1);

      u = v1;
      d = v3;
      v1 = t1;
      v3 = t3;
    } while(v3 != 0);
    v = (d - a*u)/b;
  }
  if(u_ptr != NULL) *u_ptr = u;
  if(v_ptr != NULL) *v_ptr = v;
  return d;
}

static int64_t
zz_inv_mod(
    int64_t       a,
    const int64_t p)
{
  int64_t r;
  int64_t t = ((int64_t)(a > 0) - (int64_t)(a < 0));
  a *= t;

  if(zz_gcd(a,p,&r,NULL) == 1)
  {
    r *= t;
    return (r > 0) ? r : p + r;
  }
  return 0;
}


void
pol_gen_product(
    uint16_t        *ai,
    const uint8_t   d1,
    const uint8_t   d2,
    const uint8_t   d3,
    const uint16_t  N)
{
  uint16_t rndcap = UINT16_MAX - (UINT16_MAX % N); /* assumes N is not power of 2 */
  uint16_t d = 2*((uint16_t)d1 + d2 + d3);
  uint16_t r;

  while(d > 0)
  {
    rng_uint16(&r);
    if(r < rndcap)
    {
      *ai++ = r % N;
      --d;
    }
  }

  return;
}


/* Uniform random element of pZ^n, v, such that
 * v_i + (p-1)/2 <= (q-1)/2
 * v_i - (p-1)/2 >= -(q-1)/2
 */
void
pol_unidrnd_pZ(
    int64_t          *v,
    const uint16_t   N,
    const int64_t    q,
    const int8_t     p)
{
  uint16_t i = 0;
  uint64_t r = 0;

  int64_t range = q/p;
  int64_t center = q/(2*p);

  uint64_t rndcap = (UINT64_MAX - (UINT64_MAX % range));

  while(i < N) {
    rng_uint64(&r);
    if(r < rndcap) {
      v[i] = ((int64_t)(r % range) - center) * p;
      ++i;
    }
  }

  return;
}


/*
 * Computes the inverse of a in (Z/pZ)[x]/(x^N - 1) for prime p.
 *
 * @see www.securityinnovation.com/uploads/Crypto/NTRUTech014.pdf
 *      page 3, 'Inversion in (Z/pZ)[x]/(x^N - 1)'
 */
int
pol_inv_modp(
    int64_t         *ainvp,
    const int64_t   *a,
    const uint16_t  N,
    const int64_t   p)
{
  uint16_t i;
  uint16_t k;
  uint16_t m;

  int64_t u;

  uint16_t degf;
  uint16_t degg;
  uint16_t degc;
  uint16_t degb;
  uint16_t degtmp;

  /* TODO: minimize memory usage */
  size_t scratch_len = 4*(N + 1)*sizeof(int8_t);
  int8_t *scratch = malloc(scratch_len);

  int8_t *f = scratch;
  int8_t *g = f + (N+1);
  int8_t *b = g + (N+1);
  int8_t *c = b + (N+1);
  int8_t *poltmp;

  if(!scratch)
  {
    return PQNTRU_ERROR;
  }
  memset(scratch, 0, scratch_len);

  /* f = a (mod p)*/
  degf = 0;
  for(i=0; i<N; i++) {
    f[i] = a[i] % p;
    if(f[i]) degf = i;
  }

  /* g = x^N - 1 */
  g[0] = p-1;
  g[N] = 1;
  degg = N;

  /* b(X) = 1 */
  b[0] = 1;
  degb = 0;

  /* c(X) = 0 */
  degc = 0;

  k = 0;
  while (1)
  {
    /* find smallest m such that f[m] is nonzero */
    for (m = 0; (m <= degf) && (f[m] == 0); ++m);
    if (m > degf)
    {
      free(scratch);
      return PQNTRU_ERROR; /* not invertible */
    }

    /* divide f by x^m, mul c by x^m */
    if(m > 0) {
      f = f + m;
      degf -= m;
      degc += m;
      for (i = degc; i >= m; i--)
      {
          c[i] = c[i-m];
      }
      for (i = 0; i < m; i++)
      {
          c[i] = 0;
      }
      k += m;
    }

    if(degf == 0)
    {
      break;
    }

    if(degf < degg) {
      /* Swap f and g, b and c */
      poltmp = f; f = g; g = poltmp;
      poltmp = c; c = b; b = poltmp;
      degtmp = degf; degf = degg; degg = degtmp;
      degtmp = degc; degc = degb; degb = degtmp;
    }

    u = (f[0] * zz_inv_mod(g[0], p)) % p;

    for(i = 0; i <= degg; i++)
    {
      f[i] = (f[i] - u*g[i]) % p;
    }

    if(degg == degf)
    {
      while(degf > 0 && f[degf] == 0)
      {
        --degf;
      }
    }

    for(i = 0; i <= degc; i++)
    {
      b[i] = (b[i] - u*c[i]) % p;
    }

    if(degc >= degb)
    {
      degb = degc;
      while(degb > 0 && b[degb] == 0)
      {
        --degb;
      }
    }
  }

  if(k >= N)
  {
    k -= N;
  }

  /* a^{-1} = 1/f * x^{-k} * b */
  u = zz_inv_mod(f[0], p);
  for (m=0, i = k; i < N; m++, i++)
  {
    ainvp[m] = cmod(u * b[i],  p);
  }
  for (i = 0; i < k; m++, i++)
  {
    ainvp[m] = cmod(u * b[i],  p);
  }

  shred(scratch, scratch_len);
  free(scratch);

  return PQNTRU_OK;
}


int
pol_inv_mod2(
    int64_t        *a_inv,
    const int64_t  *a,
    const uint16_t N)
{
  uint16_t i;
  uint16_t k;
  uint16_t m;

  uint16_t degf;
  uint16_t degg;
  uint16_t degc;
  uint16_t degb;
  uint16_t degtmp;

  /* TODO: minimize memory usage */
  uint16_t scratch_len = 4*(N+1);
  uint8_t *scratch = malloc(scratch_len);

  uint8_t *f = scratch;
  uint8_t *g = f + (N+1);
  uint8_t *b = g + (N+1);
  uint8_t *c = b + (N+1);
  uint8_t *poltmp;

  if(!scratch)
  {
    return PQNTRU_ERROR;
  }
  memset(scratch, 0, scratch_len);

  /* f = a (mod 2)*/
  degf = 0;
  for(i=0; i<N; i++) {
    f[i] = (uint8_t) (a[i] & 1);
    if(f[i]) degf = i;
  }

  /* g = x^N - 1 */
  g[0] = 1;
  g[N] = 1;
  degg = N;

  /* b(X) = 1 */
  b[0] = 1;
  degb = 0;

  /* c(X) = 0 */
  degc = 0;

  k = 0;

  while (1)
  {
    /* find smallest m such that f[m] is nonzero */
    for (m = 0; (m <= degf) && (f[m] == 0); ++m);
    if (m > degf)
    {
      free(scratch);
      return PQNTRU_ERROR;
    }
    if(m > 0) {
      f = f + m;
      degf -= m;
      degc += m;
      for (i = degc; i >= m; i--)
      {
          c[i] = c[i-m];
      }
      for (i = 0; i < m; i++)
      {
          c[i] = 0;
      }
      k += m;
    }

    /* if f(X) = 1, done */

    if (degf == 0)
    {
      break;
    }

    if(degf < degg) {
      /* Swap f and g, b and c */
      poltmp = f; f = g; g = poltmp;
      poltmp = c; c = b; b = poltmp;
      degtmp = degf; degf = degg; degg = degtmp;
      degtmp = degc; degc = degb; degb = degtmp;
    }

    /* f(X) += g(X)
     * might change degree of f if degg >= degf
     */

    for (i = 0; i <= degg; i++)
    {
      f[i] ^= g[i];
    }

    if(degg == degf)
    {
      while(degf > 0 && f[degf] == 0)
      {
        --degf;
      }
    }

    /* b(X) += c(X) */
    for (i = 0; i <= degc; i++)
    {
      b[i] ^= c[i];
    }

    if (degc >= degb)
    {
      degb = degc;
      while(degb > 0 && b[degb] == 0)
      {
        --degb;
      }
    }
  }

  /* a^-1 in (Z/2Z)[X]/(X^N - 1) = b(X) shifted left k coefficients */

  if (k >= N)
  {
    k = k - N;
  }

  m = 0;
  for (i = k; i < N; i++)
  {
    a_inv[m++] = (int64_t)(b[i]);
  }

  for (i = 0; i < k; i++)
  {
    a_inv[m++] = (int64_t)(b[i]);
  }

  shred(scratch, scratch_len);
  free(scratch);

  return PQNTRU_OK;
}


static void
pol_mul_indices(
    int64_t         *c,
    const int64_t   *a,
    const uint16_t  bi_P1_len,
    const uint16_t  bi_M1_len,
    const uint16_t  *bi,
    const uint16_t  N,
    int64_t         *t)
{
  //uint64_t mod_q_mask = q - 1;
  uint16_t i, j, k;

  /* t[(i+k)%N] = sum i=0 through N-1 of a[i], for b[k] = -1 */
  memset(t, 0, N * sizeof(uint64_t));

  for (j = bi_P1_len; j < bi_P1_len + bi_M1_len; j++)
  {
    k = bi[j];

    for (i = 0; k < N; ++i, ++k)
    {
      t[k] = t[k] + a[i];
    }

    for (k = 0; i < N; ++i, ++k)
    {
      t[k] = t[k] + a[i];
    }
  }

  /* t[(i+k)%N] = -(sum i=0 through N-1 of a[i] for b[k] = -1) */

  for (k = 0; k < N; k++)
  {
    t[k] = -t[k];
  }

  /* t[(i+k)%N] += sum i=0 through N-1 of a[i] for b[k] = +1 */

  for (j = 0; j < bi_P1_len; j++)
  {
    k = bi[j];

    for (i = 0; k < N; ++i, ++k)
    {
      t[k] = t[k] + a[i];
    }

    for (k = 0; i < N; ++i, ++k)
    {
      t[k] = t[k] + a[i];
    }
  }

  if(t != c)
  {
    memcpy(c, t, N*sizeof(int64_t));
  }

  return;
}


void
tri_pol_mul_product(
    int64_t         *c,
    const int64_t   *a,
    const uint16_t  b1i_len,
    const uint16_t  b2i_len,
    const uint16_t  b3i_len,
    const uint16_t  *bi,
    const uint16_t  N,
    int64_t         *t)
{
  if (b1i_len>15 || b2i_len>15 || b3i_len>15)
  {
      /* cannot use tri_pol_mul_indices when
       * maximum c_i can be greater than 31.
       */
      pol_mul_product(c, a, b1i_len,b2i_len,
                      b3i_len, bi, N, t);
      return ;
  }
  uint16_t   i;
  int64_t    *t0;
  int64_t    *t1;
  int64_t    *buf;
  t0  = t;
  t1  = t + N;
  buf = t1+ N;
  /* t2 = a * b1 */

  tri_pol_mul_indices(t0, a, b1i_len, b1i_len, bi, N, t1);

  /* t2 = (a * b1) * b2 */


  short_pol_mul_indices(t0, t0, b2i_len, b2i_len,
                  bi + (b1i_len << 1), N, buf);

  /* t = a * b3 */
  tri_pol_mul_indices(t1, a, b3i_len, b3i_len,
                  bi + ((b1i_len + b2i_len) << 1), N, buf);

  /* c = (a * b1 * b2) + (a * b3) */

  for (i = 0; i < N; i++)
  {
    c[i] = (t0[i] + t1[i]);
  }


  return;
}


void
pol_mul_product(
    int64_t         *c,
    const int64_t   *a,
    const uint16_t  b1i_len,
    const uint16_t  b2i_len,
    const uint16_t  b3i_len,
    const uint16_t  *bi,
    const uint16_t  N,
    int64_t         *t)
{
  int64_t    *t2 = t + N;
  uint16_t   i;

  /* t2 = a * b1 */

  pol_mul_indices(t2, a, b1i_len, b1i_len,
                  bi, N, t);

  /* t2 = (a * b1) * b2 */

  pol_mul_indices(t2, t2, b2i_len, b2i_len,
                  bi + (b1i_len << 1), N, t);

  /* t = a * b3 */

  pol_mul_indices(t, a, b3i_len, b3i_len,
                  bi + ((b1i_len + b2i_len) << 1), N, t);

  /* c = (a * b1 * b2) + (a * b3) */

  for (i = 0; i < N; i++)
  {
    c[i] = (t2[i] + t[i]);
  }

  return;
}

/* Space efficient Karatsuba multiplication.
 * See: Thomé, "Karatsuba multiplication with temporary space of size \le n"
 * http://www.loria.fr/~thome/files/kara.pdf
 *
 * Note: Input length should factor into b * 2^k, b <= 38
 */
static void
karatsuba(
    int64_t        *res1,   /* out - a * b in Z[x], must be length 2k */
    int64_t        *tmp1,   /*  in - k coefficients of scratch space */
    int64_t const  *a,     /*  in - polynomial */
    int64_t const  *b,     /*  in - polynomial */
    uint16_t const  k)     /*  in - number of coefficients in a and b */
{
  uint16_t i;
  uint16_t j;

  /* Grade school multiplication for small / odd inputs */
  if(k <= 38 || (k & 1) != 0)
  {
    for(j=0; j<k; j++)
    {
      res1[j] = a[0]*b[j];
    }
    for(i=1; i<k; i++)
    {
      res1[i+k-1] = 0;
      for(j=0; j<k; j++)
      {
        res1[i+j] += a[i]*b[j];
      }
    }
    res1[2*k-1] = 0;

    return;
  }

  uint16_t const p = k>>1;

  int64_t *res2 = res1+p;
  int64_t *res3 = res1+k;
  int64_t *res4 = res1+k+p;
  int64_t *tmp2 = tmp1+p;
  int64_t const *a2 = a+p;
  int64_t const *b2 = b+p;

  for(i=0; i<p; i++)
  {
    res1[i] = a[i] - a2[i];
    res2[i] = b2[i] - b[i];
  }

  karatsuba(tmp1, res3, res1, res2, p);

  karatsuba(res3, res1, a2, b2, p);

  for(i=0; i<p; i++)
  {
    tmp1[i] += res3[i];
  }

  for(i=0; i<p; i++)
  {
    res2[i]  = tmp1[i];
    tmp2[i] += res4[i];
    res3[i] += tmp2[i];
  }

  karatsuba(tmp1, res1, a, b, p);

  for(i=0; i<p; i++)
  {
    res1[i]  = tmp1[i];
    res2[i] += tmp1[i] + tmp2[i];
    res3[i] += tmp2[i];
  }

  return;
}


void
pol_mul_coefficients(
     int64_t         *c,       /* out - address for polynomial c */
     const int64_t   *a,       /*  in - pointer to polynomial a */
     const int64_t   *b,       /*  in - pointer to polynomial b */
     const uint16_t  N,        /*  in - ring degree */
     const uint16_t  padN,     /*  in - padded polynomial degree */
     const int64_t   q,        /*  in - large modulus */
     int64_t         *tmp)
{
  uint16_t i;
  int64_t *res = tmp;
  int64_t *scratch = res + 2*padN;

  karatsuba(res, scratch, a, b, padN);

  for(i=0; i<N; i++)
  {
    c[i] = cmod(res[i] + res[i+N], q);
  }
}

/* Center 'a' modulo p (an odd prime).
 * (a_i -> [-(p-1)/2, (p-1)/2]
 */
int64_t
cmod(int64_t a, int64_t p)
{
  if (a >= 0)
  {
    a %= p;
  }
  else
  {
    a = p + (a % p);
  }
  if (a > ((p-1)/2))
  {
    a -= p;
  }

  return a;
}



/* maps index = a|b into a%3|b%3 where a and b are 4 bits each */
uint8_t mod3map[256] = {
        0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,
       16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16,
       32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32,
        0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,
       16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16,
       32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32,
        0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,
       16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16,
       32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32,
        0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,
       16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16,
       32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32,
        0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,
       16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16, 17, 18, 16,
       32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32, 33, 34, 32,
        0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0};



/*
 * Both a and b are trinary; parse a into a half
 * byte polynomial; each coefficient is 0b00XX,
 * allows for 6 error-free additions; round the
 * error by mod 3 (via mod3map); 1 uint64_t addition
 * performs 16 coefficients additions.
 * */


int
pol_mul_mod_p(
     int64_t         *c,       /* out - address for polynomial c */
     const int64_t   *a,       /*  in - pointer to polynomial a */
     const int64_t   *b,       /*  in - pointer to polynomial b */
     const uint16_t  N)        /*  in - ring degree */
{
  uint16_t  i,j;
  uint64_t  *tmp;                     /* working buffer */
  uint16_t  poly_len  = (N/16+1)*8;   /* FOUR_BITS_COEFF_POLY_BYTES(P) */
  uint16_t  buf_len   = 4*poly_len;   /* store four polynomials */
  uint16_t  poly_len64= (N/16+1);     /* number of uint64_t in a poly */

  if(!( tmp= malloc(buf_len)))
  {
    return PQNTRU_ERROR;
  }
  memset(tmp, 0, buf_len);

  /* base_even = a0a1 | a2a3 | a4a5 | ... */
  uint64_t  *base_even64= (uint64_t*) tmp;
  uint8_t   *base_even  = (uint8_t*)  base_even64;

  /* base_odd  = a_N-1a0 | a1a2 | a3a4 | ... */
  uint64_t  *base_odd64 = base_even64 + poly_len64;
  uint8_t   *base_odd   = (uint8_t*)  base_odd64;

  /* accumulate results when b_i is 1 */
  uint64_t  *pos_buf64  = base_odd64 + poly_len64;
  uint8_t   *pos_buf    = (uint8_t*)  pos_buf64;
  uint16_t  pos_cnt     = 0;

  /* accumulate results when b_i is -1 */
  uint64_t  *neg_buf64  = pos_buf64 + poly_len64;
  uint8_t   *neg_buf    = (uint8_t*)  neg_buf64;
  uint16_t  neg_cnt     = 0;

  /* convert the polynomial into 4 bits poly */
  parse_4_bits_pol(base_odd, base_even,a, N);

  /* start additions; proc 2 bits from b each loop
   * performs 6 additions without cause errors; remove
   * the noise after each 6 additions.
   * */
  for (i=0;i<N/2;i++)
  {

      if (b[2*i]==1)
      {
          pos_cnt++;
          for (j=0;j<poly_len64;j++)
          {
              pos_buf64[j] += base_even64[j];
          }
          if (pos_cnt==6)
          {
              pos_cnt=0;
              for (j=0;j<poly_len;j++)
                  pos_buf[j]  = mod3map[pos_buf[j]];
         }
      }
      else if (((int)b[2*i])==-1 || b[2*i]==2)
      {
          neg_cnt++;
          for (j=0;j<poly_len64;j++)
          {
              neg_buf64[j] += base_even64[j];
          }
          if (neg_cnt==6)
          {
              neg_cnt=0;
              for (j=0;j<poly_len;j++)
                  neg_buf[j]  = mod3map[neg_buf[j]];
          }
      }
        if (b[2*i+1]==1)
        {
            pos_cnt++;
            for (j=0;j<poly_len64;j++)
            {
                pos_buf64[j] += base_odd64[j];
            }
            if (pos_cnt==6)
            {
                pos_cnt=0;
                for (j=0;j<poly_len;j++)
                    pos_buf[j]  = mod3map[pos_buf[j]];
            }
        }
        else if (((int)b[2*i+1])==-1 || b[2*i+1]==2)
        {
            neg_cnt++;
            for (j=0;j<poly_len64;j++)
            {
                neg_buf64[j] += base_odd64[j];
            }
            if (neg_cnt==6)
            {
                neg_cnt=0;
                for (j=0;j<poly_len;j++)
                    neg_buf[j]  = mod3map[neg_buf[j]];
            }
        }

      /* shift both polynomials */
      double_shift_poly(base_even,N);
      double_shift_poly(base_odd, N);

  }

  /* taking care of the last coefficient of b*/
  if (b[N-1]==1)
  {
      pos_cnt++;
      for (j=0;j<poly_len64;j++)
      {
          pos_buf64[j] += base_even64[j];
      }
      if (pos_cnt==6)
      {
          pos_cnt=0;
          for (j=0;j<poly_len;j++)
              pos_buf[j]  = mod3map[pos_buf[j]];
      }
  }
  else if (((int)b[N-1])==-1 || b[N-1]==2)
  {
      neg_cnt++;
      for (j=0;j<poly_len64;j++)
      {
          neg_buf64[j] += base_even64[j];
      }
      if (neg_cnt==6)
      {
          neg_cnt=0;
          for (j=0;j<poly_len;j++)
              neg_buf[j]  = mod3map[neg_buf[j]];
      }
  }

  /* put final result into c */
  memset(c,0,N*8);

  for( i=0;i<N/2;i++)
  {
      c[2*i]    = cmod((pos_buf[i]& 0b1111)
                      -(neg_buf[i]& 0b1111),3);
      c[2*i+1]  = cmod(((pos_buf[i]& 0b11110000)>>4)
                      -((neg_buf[i]& 0b11110000)>>4),3);
  }
  c[N-1] = cmod((pos_buf[N/2]& 0b1111)
               -(neg_buf[N/2]& 0b1111),3);


  free(tmp);
  return 0;
}

/*
 * Functions for two bytes trinary polynomials
 */
/*
 * a is short |a|<30; b is index (product form)
 * compute c = a* b.coeff.
 * convert a *int64_t poly into *int16_t poly
 * where -1 = 0b0000 0111 1111 1111, this allows
 * for 31 error free additions; no need to handle
 * the error as we do maxmimum df3*2 = 30 additions;
 * use uint64_t additions - 1 uint64_t addition
 * adds 5 coefficients; the coefficient of
 * result c should be in between (-1023, 1023);
 * */

void
short_pol_mul_indices(
    int64_t         *c,
    const int64_t   *a,
    const uint16_t  bi_len,
    const uint16_t  bi_M1_len,
    const uint16_t  *bi,
    const uint16_t  N,
    uint64_t         *t)
{
    uint16_t    i,j;
    int8_t      t1,t2;
    uint16_t    poly_len64;
    uint16_t    *base_poly;
    uint64_t    *base_poly64;
    uint16_t    *rot_poly;
    uint64_t    *rot_poly64;
    uint16_t    *pos_res_poly;
    uint64_t    *pos_res_poly64;
    uint16_t    *neg_res_poly;
    uint64_t    *neg_res_poly64;
    uint16_t    *tmp;
    tmp = t;
    poly_len64  =   N/4+1;      /*  number of uint64_t in a poly */

    memset(tmp,  0, 4*poly_len64*sizeof(uint64_t));

    base_poly64     = (uint64_t*) tmp;
    base_poly       = (uint16_t*) base_poly64;

    rot_poly64      = base_poly64 + poly_len64;
    rot_poly        = (uint16_t*) rot_poly64;

    pos_res_poly64  = rot_poly64  + poly_len64;
    pos_res_poly    = (uint16_t*) pos_res_poly64;

    neg_res_poly64  = pos_res_poly64  + poly_len64;
    neg_res_poly    = (uint16_t*) neg_res_poly64;

    parse_two_bytes_pol(base_poly,a,N);

    /* handling positive coefficients */
    for (i=0;i<bi_len;i++)
    {
        shift_two_bytes_poly_k(rot_poly,base_poly, bi[i], N);

        for (j=0;j<poly_len64;j++)
        {
            pos_res_poly64[j] += rot_poly64[j];
        }
    }

    /* handling negative coefficients */
    for (i=bi_len;i<2*bi_len;i++)
    {
        shift_two_bytes_poly_k(rot_poly,base_poly, bi[i], N);

        for (j=0;j<poly_len64;j++)
        {
            neg_res_poly64[j] += rot_poly64[j];
        }
    }

    /* subtract neg from pos */

    for (j=0;j<N;j++)
    {
        t1  =  pos_res_poly[j] & 0x07ff;
        if (t1>1024)
            t1 -=2048;
        else if (t1<-1024)
            t1 +=2048;

        t2  =  neg_res_poly[j] & 0x07ff;
        if (t2>1024)
            t2 -=2048;
        else if (t2<-1024)
            t2 +=2048;

        c[j] = t1-t2;
    }
    return ;
}

/*
 * convert a *int64_t poly into *int16_t poly
 * where |a| < 30; using 11 bits to represent
 * a and 5 bits of 0 to guard.
 * */
void
parse_two_bytes_pol(
    uint16_t        *out,
    const int64_t   *in,
    const uint16_t  N)
{
    uint16_t i;
    uint16_t  *buf;
    buf = out;
    memset(buf, 0, N);
    for (i=0;i<N;i++)
    {
        buf[i] = in[i]&0x07FF;
    }
    return ;
}

/*
 * each coefficient of poly is two bytes;
 * shift poly by 1 byte/coefficient
 */
void
shift_two_bytes_poly(
    uint16_t         *p,
    const uint16_t  N)
{
    uint16_t tmp;
    tmp =   p[N-1];
    memcpy(p+1, p, 2*(N-1));
    p[0] = tmp;
}
/*
 * each coefficient of poly is two bytes;
 * shift poly by "index" coefficients
 * */
void
shift_two_bytes_poly_k(
    uint16_t        *out,
    const uint16_t  *in,
    const uint16_t  index,
    const uint16_t  N)
{
    memcpy(out,         in+N-index, 2*index);
    memcpy(out+index,   in,         2*(N-index));
    return;
}

/*
 * Functions for full byte trinary polynomials
 */

/*
 * each coefficient of poly is a byte;
 * shift poly by 1 byte/coefficient
 */
void
shift_poly(
    uint8_t         *p,
    const uint16_t  N)
{
    uint8_t tmp;
    tmp =   p[N-1];
    memcpy(p+1, p, N-1);
    p[0] = tmp;
}

/*
 * each coefficient of poly is a byte;
 * shift poly by "index" byte/coefficient
 * */
void
shift_poly_k(
    uint8_t         *out,
    const uint8_t   *in,
    const uint16_t  index,
    const uint16_t  N)
{
    memcpy(out,         in+N-index, index);
    memcpy(out+index,   in,         N-index);
    return;
}


/*
 * a is trinary; b is index (product form)
 * compute c = a* b.coeff.
 * convert a *int64_t poly into *int8_t poly
 * where -1 = 0b00111111, this allows for 3
 * error free additions; round off the error
 * by &0b00111111 (&x3F) after three additions;
 * use uint64_t additions - 1 uint64_t addition
 * adds 8 coefficients; the coefficient of
 * result c should be in between (-31, 31);
 * */

void
tri_pol_mul_indices(
    int64_t         *c,
    const int64_t   *a,
    const uint16_t  bi_len,
    const uint16_t  bi_M1_len,
    const uint16_t  *bi,
    const uint16_t  N,
    uint64_t         *t)
{

    uint16_t    i,j;
    int8_t      t1,t2;
    uint16_t    poly_len64;
    uint8_t     *base_poly;
    uint64_t    *base_poly64;
    uint8_t     *rot_poly;
    uint64_t    *rot_poly64;
    uint8_t     *pos_res_poly;
    uint64_t    *pos_res_poly64;
    uint8_t     *neg_res_poly;
    uint64_t    *neg_res_poly64;
    uint8_t     *tmp;
    tmp = t;
    poly_len64  =   N/8+1;      //  number of uint64_t in a poly

    memset(tmp,  0, 4*poly_len64*sizeof(uint64_t));

    base_poly64     = (uint64_t*) tmp;
    base_poly       = (uint8_t*)  base_poly64;

    rot_poly64      = base_poly64 + poly_len64;
    rot_poly        = (uint8_t*)  rot_poly64;

    pos_res_poly64  = rot_poly64  + poly_len64;
    pos_res_poly    = (uint8_t*)  pos_res_poly64;

    neg_res_poly64  = pos_res_poly64  + poly_len64;
    neg_res_poly    = (uint8_t*)  neg_res_poly64;

    parse_pol(base_poly,a,N);

    /* handling positive coefficients */
    for (i=0;i<bi_len;i++)
    {
        shift_poly_k(rot_poly,base_poly, bi[i], N);

        for (j=0;j<poly_len64;j++)
        {
            pos_res_poly64[j] += rot_poly64[j];
        }
        if (i%3==2)
        {
            for (j=0;j<poly_len64;j++)
            {
                pos_res_poly64[j] &= 0x3f3f3f3f3f3f3f3f;
            }
        }
    }

    /* handling negative coefficients */
    for (i=bi_len;i<2*bi_len;i++)
    {
        shift_poly_k(rot_poly,base_poly, bi[i], N);

        for (j=0;j<poly_len64;j++)
        {
            neg_res_poly64[j] += rot_poly64[j];
        }

        if (i%3==2)
        {
            for (j=0;j<poly_len64;j++)
            {
                neg_res_poly64[j] &= 0x3f3f3f3f3f3f3f3f;
            }
        }
    }

    /* subtract neg from pos */

    for (j=0;j<N;j++)
    {
        t1  =  pos_res_poly[j] & 0x3f;
        if (t1>32)
            t1 -=64;
        else if (t1<-32)
            t1 +=64;

        t2  =  neg_res_poly[j] & 0x3f;
        if (t2>32)
            t2 -=64;
        else if (t2<-32)
            t2 +=64;

        c[j] = t1-t2;
    }
    return ;
}

/*
 * convert a *int64_t poly into *int8_t poly
 * where -1 = 0b00111111
 * */
void
parse_pol(
    uint8_t         *out,
    const int64_t   *in,
    const uint16_t  N)
{
    uint16_t i;
    uint8_t  *buf;
    buf = out;
    memset(buf, 0, N);
    for (i=0;i<N;i++)
    {
        switch (in[i])
         {
             case(1):    buf[i] = 0b00000001; break;
         /*  case(0):    buf[i] = 0b00000000; break; */
             case(-1):   buf[i] = 0b00111111; break;
         }
    }
    return ;
}


/*
 * Functions for half byte trinary polynomials
 */
void
double_shift_poly(
    uint8_t         *p,
    const uint16_t  N)
{
    uint8_t tmp1;
    uint8_t tmp2;
    uint16_t len = N/2+1;
    tmp1 =   p[len-1];          // a_N-1,0
    tmp2 =   p[len-2];          // a_N-3,a_N-2

    memmove(p+1, p, len-1);     // shift memory by 1 byte

    p[0]                        // reset the first byte
                                // should be a_N-2|a_N-1
         = ((tmp1 & 0b1111)     // last 4 bits of tmp1 (a_N-1)
            <<4)                //   shift left
         + ((tmp2 & 0b11110000) // top 4 bits of tmp2 (a_N-2)
            >>4) ;              //   shift right
    return ;
}


void
parse_4_bits_pol(
    uint8_t         *out_odd,
    uint8_t         *out_even,
    const int64_t   *in,
    const uint16_t  N)
{
    uint16_t i,len;
    uint8_t  *buf_odd;
    uint8_t  *buf_even;
    buf_odd     = out_odd;
    buf_even    = out_even;
    len         = N/2+1;
    memset(buf_odd, 0, len);
    memset(buf_even, 0, len);

    for (i=0;i<N/2;i++)
    {
    /* allocate even array a0a1 | a2a3 | a4a5 | ... | aN-3aN-2 | EMPTY */

        if (in[2*i]==1&&in[2*i+1]==1)
            buf_even[i] = 0b00010001;
        if (in[2*i]==1&&in[2*i+1]==0)
            buf_even[i] = 0b00000001;
        if (in[2*i]==1&&in[2*i+1]==-1)
            buf_even[i] = 0b00100001;

        if (in[2*i]==0&&in[2*i+1]==1)
            buf_even[i] = 0b00010000;
     /* if (in[2*i]==0&&in[2*i+1]==0)
            buf_even[i] = 0b00000000; */
        if (in[2*i]==0&&in[2*i+1]==-1)
            buf_even[i] = 0b00100000;

        if (in[2*i]==-1&&in[2*i+1]==1)
            buf_even[i] = 0b00010010;
        if (in[2*i]==-1&&in[2*i+1]==0)
            buf_even[i] = 0b00000010;
        if (in[2*i]==-1&&in[2*i+1]==-1)
            buf_even[i] = 0b00100010;

    /* allocate odd array EMPTY | a1a2 | a3a4 | a5a6 | ... | aN-2aN-1 */

        if (in[2*i+1]==1&&in[2*i+2]==1)
            buf_odd[i+1] = 0b00010001;
        if (in[2*i+1]==1&&in[2*i+2]==0)
            buf_odd[i+1] = 0b00000001;
        if (in[2*i+1]==1&&in[2*i+2]==-1)
            buf_odd[i+1] = 0b00100001;

        if (in[2*i+1]==0&&in[2*i+2]==1)
            buf_odd[i+1] = 0b00010000;
    /*  if (in[2*i+1]==0&&in[2*i+2]==0)
            buf_odd[i+1] = 0b00000000;  */
        if (in[2*i+1]==0&&in[2*i+2]==-1)
            buf_odd[i+1] = 0b00100000;

        if (in[2*i+1]==-1&&in[2*i+2]==1)
            buf_odd[i+1] = 0b00010010;
        if (in[2*i+1]==-1&&in[2*i+2]==0)
            buf_odd[i+1] = 0b00000010;
        if (in[2*i+1]==-1&&in[2*i+2]==-1)
            buf_odd[i+1] = 0b00100010;

    }
    /* end pad the even array with aN-1,0 */
    if (in[N-1] ==1)
        buf_even[N/2] = 0b00000001;
    /* if (in[N-1] ==0)
        buf_even[N/2] = 0b00000000; */
    if (in[N-1] ==-1)
        buf_even[N/2] = 0b00000010;

    /* front pad the odd array with aN-1,a0 */
    if (in[0]==1&&in[N-1]==1)
        buf_odd[0] = 0b00010001;
    if (in[0]==1&&in[N-1]==0)
        buf_odd[0] = 0b00010000;
    if (in[0]==1&&in[N-1]==-1)
        buf_odd[0] = 0b00010010;

    if (in[0]==0&&in[N-1]==1)
        buf_odd[0] = 0b00000001;
/*  if (in[0]==0&&in[N-1]==0)
        buf_odd[i+1] = 0b00000000;  */
    if (in[0]==0&&in[N-1]==-1)
        buf_odd[0] = 0b00000010;

    if (in[0]==-1&&in[N-1]==1)
        buf_odd[0] = 0b00100001;
    if (in[0]==-1&&in[N-1]==0)
        buf_odd[0] = 0b00100000;
    if (in[0]==-1&&in[N-1]==-1)
        buf_odd[0] = 0b00100010;

    return;
}

