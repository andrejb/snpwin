/*
 * read_data.c - routine for reading from source SNP file
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
#include <errno.h>
#include "read_data.h"

/*
 * This routine reads from a preformatted file with SNP zygosity information
 * (1 - homozygote, 2 - heterozygote), and stores in the **data pointer. Also,
 * direct pointers for each SNP data are stored in the ***snpxind vector.
 *
 * For more information, see the ../README file.
 */
void
read_data (char ***snpxind, char **data, FILE * stream, int snp_num,
	   int ind_num)
{
  char buf[STRLEN_MAX];
  int pos, p, s;
  char val;

  /* memory space for snp pointers */
  if ((*snpxind = (char **) malloc (snp_num * sizeof (char *))) == NULL)
    {
      fprintf (stderr, "Error allocating memory for data.\n");
      exit (2);
    }

  /* memory space for actual data */
  if ((*data = (char *) malloc (snp_num * ind_num * sizeof (char))) == NULL)
    {
      fprintf (stderr, "Error allocating memory for data.\n");
      exit (2);
    }

  /* links snpxind to pdata entries */
  p = 0;
  for (s = 0; s < snp_num; s++)
    {
      (*snpxind)[s] = (*data + p);
      p += ind_num;
    }

  s = 0;
  while ((fgets (buf, STRLEN_MAX, stream) != NULL) &&
	 buf[0] != '\0' && s < snp_num)
    {
      /*printf("read a string: %s\n", buf); */
      pos = 0;

      /* advances to third <TAB> */
      while (buf[pos] != '\t')
	pos++;
      pos++;
      while (buf[pos] != '\t')
	pos++;
      pos++;
      while (buf[pos] != '\t')
	pos++;


      /* reads ind_num SNPs */
      for (p = 0; p < ind_num; p++)
	{
	  /* this ignores unwanted data (such as ^M and 'A's from 'NA's) */
	  while (buf[pos] != '0' && buf[pos] != '1'
		 && buf[pos] != '2' && buf[pos] != 'N')
	    pos++;

	  /* this stores read values, '3' means NA. */
	  val = buf[pos];
	  if (val == '0' || val == '1' || val == '2')
	    ((*snpxind)[s])[p] = val;
	  else
	    {
	      ((*snpxind)[s])[p] = '3';
	    }
	  pos++;
	}
      s++;
    }

}

/*
 * routine for openning a file for writing.
 */
FILE *
open_file (char *fname, char *mode)
{
  FILE *fp;
  if ((fp = fopen (fname, mode)) == NULL)
    {
      fprintf (stderr, "Error (errno: %d) opening file '%s' for writing.",
	       errno, fname);
      exit (1);
    }
  return fp;
}
