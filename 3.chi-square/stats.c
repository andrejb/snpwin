/*
 * stats.c - calculate min, max and mean windows sizes (with stddev)
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

#include "../common/util.h"

#define STRLEN_MAX 30

#define max(x,y) (x > y ? x : y)
#define min(x,y) (x < y ? x : y)


/*
 * routine for opening a file for writing.
 */
FILE *
open_file (char *fname)
{
  FILE *fp;
  if ((fp = fopen (fname, "w")) == NULL)
    {
      fprintf (stderr, "Error (errno: %d) opening file '%s' for writing.",
	       errno, fname);
      exit (1);
    }
  return fp;
}


/*
 * calculate statistics about windows (min, max and mean size, and standard
 * deviation.
 */
int
main (int argc, char *argv[])
{

  int win, start, end, max, len, i, min;
  float mean, sum;
  char buf[STRLEN_MAX + 1];
  int sizes[48697];


  /* read stdin */
  win = 0;
  max = 0;
  sum = 0;
  min = 10;
  while ((fgets (buf, STRLEN_MAX, stdin) != NULL) && buf[0] != '\0')
    {
      sscanf (buf, "%d %d", &start, &end);
      len = end - start + 1;
      max = max (max, len);
      min = min (min, len);
      sizes[win] = len;
      sum += len;
      win++;
    }

  mean = sum / win;

  printf ("Sum:  %f\n", ((float) sum));
  printf ("Mean: %f\n", mean);
  printf ("Max:  %d\n", max);
  printf ("Min:  %d\n", min);

  sum = 0;
  for (i = 0; i < win; i++)
    {
      sum += fabs (mean - sizes[i]);
    }
  printf ("Standard deviation: %f\n", sum / win);

  return 0;
}
