/*
 * util.c - library routines for SNP neighborhood determination
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "util.h"
#include "tree.h"

/*
 * limit_neighborhoods: verifies maximum neighborhood length given available
 * data.
 */
void
limit_neighborhoods (pair_t * neighborhoods,
		     char **snpxind, int snp, int ind_num, int lmax, int rmax)
{
  int l, r, i;
  for (i = 0; i < ind_num; i++)
    {
      l = 0;
      while (l < lmax && snptoi (snpxind[snp - l - 1][i]) != -1)
	l++;
      r = 0;
      while (r < rmax && snptoi (snpxind[snp + r + 1][i]) != -1)
	r++;
      neighborhoods[i].l = l;
      neighborhoods[i].r = r;
    }
  return;
}


/*
 * malloc_or_fail:
 */
void *
malloc_or_fail (size_t size)
{
  void *p;
  if ((p = (void *) malloc (size)) == NULL)
    {
      fprintf (stderr, "Error allocating memory (errno: %d).\n", errno);
      exit (1);
    }
  return p;
}
