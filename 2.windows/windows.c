/*
 * windows.c - parse neighborhood results and generate some info about windows
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
#include "../common/read_data.h"

/* output file names */
#define FILE_SUPERPOS       "../results/windows-superpos.txt"
#define FILE_INTERVALS      "../results/windows-intervals.txt"
#define FILE_RATIO          "../results/windows-ratio.txt"
#define FILE_WINDOWS        "../results/windows-windows.txt"
#define FILE_INTERVAL_SIZES "../results/windows-interval-sizes.txt"

#define max(x,y) (x > y ? x : y)

/*
 * the main routine allocates memory, opens files for reading and performs
 * analysis of the neighborhoods result file for calculating windows.
 */
int
main (int argc, char *argv[])
{

  /*FILE *superpos, *intervals, *ratio, *interval_sizes;*/
  FILE *windows;
  int snp, l, r, total;
  pair_t *neighborhoods;
  int *superposition;
  float iratio;
  char buf[STRLEN_MAX + 1];
  float *position, pos;
  int rightmost, j, window, leftmost;

  /* decide on total number of SNPs based on command line argument */
  if (argc < 2)
    total = SNP_NUM;
  else
    total = atoi (argv[1]);
  printf ("total: %d\n", total);

  /* allocate needed space */
  superposition = malloc_or_fail (total * sizeof (int));
  neighborhoods = malloc_or_fail (total * sizeof (pair_t));
  position = malloc_or_fail (total * sizeof (float));

  /* read neighborhood data from stdin */
  snp = 0;
  while (snp < total &&
	 (fgets (buf, STRLEN_MAX, stdin) != NULL) && buf[0] != '\0')
    {
      sscanf (buf, "%f: (%d, %d) %f", &pos, &l, &r, &iratio);
      position[snp] = pos;
      neighborhoods[snp].l = l;
      neighborhoods[snp].r = r;
      neighborhoods[snp].ratio = iratio;

      /*
         if (snp == 500179)
         {
         printf ("    %f: (%d, %d) %f\n", pos, l, r, iratio);
         printf ("    %f: (%d, %d) %f\n", position[snp],
         neighborhoods[snp].l, neighborhoods[snp].r,
         neighborhoods[snp].ratio);
         }
       */

      snp++;
    }
  printf ("Read %d lines.\n", snp);

  /* open files for writing */
  /*superpos = open_file (FILE_SUPERPOS, "w");
  intervals = open_file (FILE_INTERVALS, "w");
  interval_sizes = open_file (FILE_INTERVAL_SIZES, "w");
  ratio = open_file (FILE_RATIO, "w");*/
  windows = open_file (FILE_WINDOWS, "w");


  /* initialize superposition vector with zeros */
  for (snp = 0; snp < total; snp++)
    superposition[snp] = 0;

  /* iterate through SNPs to calculate windows */
  rightmost = 1;
  leftmost = 0;
  for (snp = 0; snp < total; snp++)
    {
      /* window determination */
      rightmost--;
      rightmost = max (rightmost, neighborhoods[snp].r);
      if (rightmost == 0)	/* got a potential window! */
	{
	  window = 1;
	  for (j = snp + 1; j < snp + 20; j++)
	    if ((j - neighborhoods[j].l) <= snp)
	      {
		window = 0;
		break;
	      }

	  /* if a window was found, write to output  */
	  if (window == 1)
	    {
	      fprintf (windows, "%.0f %.0f\n", position[leftmost],
		       position[snp]);
	      leftmost = snp + 1;
	    }
	}

      /* calculates superposition between snps */
      for (l = 1; l <= neighborhoods[snp].l; l++)
	superposition[snp - l]++;
      for (r = 1; r <= neighborhoods[snp].r; r++)
	{
	  if ((snp + r) == total)	/* take care with limits to the right */
	    break;
	  superposition[snp + r]++;
	}
    }

  /* print analysis results */
  //for (snp = 0; snp < total; snp++)
  //  {
  //    /* superpos.txt */
  //    fprintf (superpos, "%.0f %d\n", position[snp], superposition[snp]);
  //    /* intervals.txt */
  //    fprintf (intervals,
  //             "%.0f %d -%d\n",
  //             position[snp], neighborhoods[snp].l, neighborhoods[snp].r);
  //    /* ratio.txt */
  //    fprintf (ratio, "%.0f %f\n", position[snp],
  //             -2 * log (neighborhoods[snp].ratio));
  //    /* interval-sizes.txt */
  //    fprintf (interval_sizes, "%.0f %d\n", position[snp],
  //             neighborhoods[snp].l + neighborhoods[snp].r);
  //  }

  /* close file handlers */
  fclose (stdin);
  /*fclose (superpos);
  fclose (intervals);
  fclose (interval_sizes);
  fclose (ratio);*/
  fclose (windows);
  return 0;
}
