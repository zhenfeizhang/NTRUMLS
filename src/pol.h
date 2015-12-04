/*
 * CPQREF/pol.h
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

#ifndef CPQREF_POL_H_
#define CPQREF_POL_H_


void
pol_gen_product(
    uint16_t        *ai,
    const uint8_t   d1,
    const uint8_t   d2,
    const uint8_t   d3,
    const uint16_t  N);


void
pol_unidrnd_pZ(
    int64_t          *v,
    const uint16_t   N,
    const int64_t    q,
    const int8_t     p);


int
pol_inv_modp(
    int64_t         *ainvp,
    const int64_t   *a,
    const uint16_t  N,
    const int64_t   p);


int
pol_inv_mod2(
    int64_t        *a_inv,
    const int64_t  *a,
    const uint16_t N);


void
pol_mul_product(
    int64_t         *c,
    const int64_t   *a,
    const uint16_t  b1i_len,
    const uint16_t  b2i_len,
    const uint16_t  b3i_len,
    const uint16_t  *bi,
    const uint16_t  N,
    int64_t         *t);


void
pol_mul_coefficients(
     int64_t         *c,
     const int64_t   *a,
     const int64_t   *b,
     const uint16_t  N,
     const uint16_t  padN,
     const  int64_t  q,
     int64_t         *tmp);



int64_t cmod(const int64_t a, const int64_t p);


/*
 * Functions for half byte trinary polynomials
 * each coefficient is 4 bits;
 */

/*
 * shift coeffients of a poly twice; e.g.,
 * in  = a0   + a1x   + a2x^2 ... + aN-1x^N-1
 * out = aN-2 + aN-1x + a0x^2 ... + aN-3x^N-3
 */
void
double_shift_poly(
    uint8_t         *p,
    const uint16_t  N);

/*
 * convert a *int64_t poly into half byte poly
 * where -1 = 0b0010; 0 = 0b0000; 1 = 0b0001
 * outputs two polynomials:
 * odd  = a1   + a2x   + a3x^2 ... + a0x^N-1
 * even = a0   + a1x   + a2x^2 ... + aN-1x^N-1
 * */

void
parse_4_bits_pol(
    uint8_t         *out_odd,
    uint8_t         *out_even,
    const int64_t   *in,
    const uint16_t  N);

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
     const uint16_t  N);       /*  in - ring degree */


/*
 * Functions for full byte trinary polynomials
 */

/*
 * each coefficient of poly is a byte;
 * shift poly by "index" byte/coefficient
 * */
void
shift_poly_k(
    uint8_t         *out,
    const uint8_t   *in,
    const uint16_t  index,
    const uint16_t  N);


/*
 * a is trinary; b is index (product form)
 * compute c = a* b.coeff.
 * convert a *int64_t poly into *int8_t poly
 * where -1 = 0b00111111, this allows for 3
 * error free additions; round off the error
 * by &0b00111111 (&03F) after three additions;
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
    uint64_t         *t);

/*
 * convert a *int64_t poly into *int8_t poly
 * where -1 = 0b00111111
 * */
void
parse_pol(
    uint8_t         *out,
    const int64_t   *in,
    const uint16_t  N);
#endif //CPQREF_POL_H_
