/*
 * util.h - definitions for util.c
 *
 * Copyright (c) 2009 André Jucovsky Bianchi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _util_h
#define _util_h

#include "tree.h"

/* maximum neighborhood length */
/*#define MAX_NLEN 20*/

/* the penalizing constant is (|A|-1)(|A|^k)(log(n)) */
#define PENALIZING_TERM(k, n) (1 * pow(3, k) * (logtable[n][1]))

/* neighobrhood struct */
typedef struct pair
{
  int l;
  int r;
  float ratio;
} pair_t;

/* function prototypes */
void limit_neighborhoods (pair_t * neighborhoods,
			  char **snpxind, int snp, int ind_num, int lmax, int
			  rmax);
void *
malloc_or_fail (size_t size);

#endif /* _util_h */
