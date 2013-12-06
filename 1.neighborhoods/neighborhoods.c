/*
 * neighborhoods.c - estimate SNP neighborhoods
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

#include "../common/read_data.h"
#include "../common/util.h"
#include "../common/tree.h"

//#define SNP_DEBUG

/* precalculated log table to avoid having to calculate logs every time */
double **logtable = NULL;


/* function prototypes */

void estimate_neighborhood (char **snpxind, int snp_num, int ind_num, int snp,
			    pair_t * p, int max_nlen);
void calc_likelihood (tree_type T, double *likelihood, int k);
void bootstrap_memory (double ***likelihood, pair_t ** neighborhoods, int
		       max_nlen, int ind_num);

/* main: entry point for SNP neighborhood estimation.
 *
 * The main routine controls reading from input, estimating each SNP
 * neighborhood with maximum likelihood for each SNP, and writing to output.
 */
int main (int argc, char *argv[]) {

  FILE *infile;			/* pointer for input infile */
  char *data;			/* table to hold raw SNP data */
  char **snpxind;		/* indexed pointers for raw SNP data */
  int snp, start;		/* auxiliary variables */
  pair_t neighborhoods[SNP_NUM];	/* neighborhood results */
  int max_nlen;			/* maximum neighborhood length */

  /*int superposition[SNP_NUM]; */

  /* read the input infile */
  infile = open_file (argv[1], "r");
  read_data (&snpxind, &data, infile, SNP_NUM,	/* total number of SNPs */
	     IND_NUM);		/* number of individuals, i.e. sample size */

  /* get starting point */
  if (argc == 3)
    start = atoi (argv[2]);
  else
    start = 0;

  /* calculate maximum neighborhood length */
  max_nlen = ceil (log10 (SNP_NUM) / log10(3));
  fprintf (stderr, "Maximum neighborhood length: %d\n", max_nlen);

  /* determine best neighborhood for each snp ... */
  for (snp = start; snp < SNP_NUM; snp++)
    {
      estimate_neighborhood (snpxind,	/* memory space for data */
			     SNP_NUM, IND_NUM, snp, &(neighborhoods[snp]),
			     max_nlen);
      /* ... and print the results. */
      printf ("%d: (%d, %d) %f\n",
	      snp,
	      neighborhoods[snp].l,
	      neighborhoods[snp].r, neighborhoods[snp].ratio);
      if ((snp % 10) == 0)
	fflush (stdout);
    }

  fclose (infile);
  return 0;
}

/*
 * estimate_neighborhood: estimate left and right neighborhoods for a SNP.
 *
 * This routine iterates through all possible sizes for left and right
 * neighborhood, calculates the likelihood for each of them, and finally
 * stores the values for left and right neighborhoods that achieved the
 * maximum likelihood given case and control samples.
 */
void
estimate_neighborhood (char **snpxind, int snp_num, int ind_num, int snp,
		       pair_t * p, int max_nlen)
{
  int ind, l, r, lmax, rmax, skip, start;
  tree_type T;
  static double **likelihood = NULL;
  static pair_t *neighborhoods = NULL;

  /* allocate memory, calculate tables, etc. */
  bootstrap_memory (&likelihood, &neighborhoods, max_nlen, ind_num);

  /* determine left and right limits */
  lmax = max_nlen;
  if (snp - lmax < 0)
    lmax = snp;
  rmax = max_nlen;
  if (snp + rmax > snp_num - 1)
    rmax = snp_num - snp - 1;

  /* limit actual neighborhoods to be considered based on NA entries */
  limit_neighborhoods (neighborhoods, snpxind, snp, ind_num, lmax, rmax);

  /* estimation of neighborhood with maximum likelihood */
  for (l = 0; l <= lmax; l++)	/* all possible left neighborhoods */
    for (r = 0; r <= rmax; r++)	/* all possible right neighborhoods */
      {
	/*printf("Calculating for SNP %d with l=%d and r=%d\n", snp, l, r); */
	T = make_empty_tree ();
	/* look for all individuals */
	for (ind = 0; ind < ind_num; ind++)
	  {
	    /* depending on the available info */
	    if (snptoi (snpxind[snp][ind]) != -1	/* SNP `i` is not NA */
		&& l <= neighborhoods[ind].l	/* lefs side contains no NA */
		&& r <= neighborhoods[ind].r)	/* right side contains no NA */
	      {
		/* indicates if we'll have to skip a position */
		if (l == 0 || r == 0)
		  skip = -1;
		else
		  skip = l;

		/* calculate start point of insertion */
		start = snp - l - 1;
		if (l == 0)
		  start = snp;

		/* insert the string on the tree. */
		insert_string (T,	/* the tree */
			       snpxind,	/* SNP data */
			       ind,	/* which individual's SNPs are being inserted */
			       start,	/* starting point of insertion */
			       l + r,	/* total length of neighborhood candidate */
			       snptoi (snpxind[snp][ind]),	/* central SNP value */
			       skip);	/* avoid inserting the central SNP */
	      }
	  }
	/*print_tree(T, 0, 0); */

	/* calculate the likelihood of the neighborhood with `l` SNPs to the
	   left and `r` SNPs to the right. */
        likelihood[l][r] = 0;
	calc_likelihood (T, &(likelihood[l][r]), l + r);

	/* penalize likelihood */
#ifdef SNP_DEBUG
	fprintf (stderr, "SNP %d: #%d likelihood(l=%d, r=%d) = %f - %f\n",
		snp, T->counter[3],
		l, r, likelihood[l][r], PENALIZING_TERM ((l + r),
							 T->counter[3]));
	fprintf (stderr, "not penalized likelihood - l: %d r: %d v: %f\n", l, r,
		likelihood[l][r]);
#endif
	likelihood[l][r] -= PENALIZING_TERM ((l + r), T->counter[3]);
#ifdef SNP_DEBUG
	fprintf (stderr, "penalized likelihood     - l: %d r: %d v: %f\n", l, r,
		likelihood[l][r]);
#endif
	free_tree (T);
      }


  /* find `l` and `r` with maximum likelihood */
  p->l = 0;
  p->r = 0;
  for (l = 0; l <= lmax; l++)	/* all possible left neighborhoods */
    for (r = 0; r <= rmax; r++)	/* all possible right neighborhoods */
      if (likelihood[l][r] > likelihood[p->l][p->r])
	{
	  p->l = l;
	  p->r = r;
	}
  /*printf("Chosen neighborhood: l=%d, r=%d: %f\n", p->l, p->r, likelihood[p->l][p->r]); */

  /* ratio between likelihoods of chosen neighborhood and no neighborhood at all */
  p->ratio = likelihood[0][0] / likelihood[p->l][p->r];
  return;
}

/*
 * calc_likelihood: calculate the likelihood of one neighborhood for one SNP.
 *
 * This routine descends on the tree containing the SNP information for all
 * individuals (considering specific `l` and `r` candidates for left and right
 * neighborhood values centered on one specific SNP) and calculates its
 * likelihood.
 */
void
calc_likelihood (tree_type T, /* tree to analyse */
                  double *likelihood, /* likelihood data will be stored here */
                  int k) /* depth to descend before calculating */
{
  int a;
  double s;

  /* if this node is a leaf, then calculate likelihood of this choice for neighborhood */
  if (k == 0)
    for (a = 0; a <= 2; a++)
      {
	if (T->counter[a] != 0)
	  {
	    s =
	      ((double) T->counter[a]) *
	      (logtable[T->counter[a]][T->counter[3]]);
/*#ifdef SNP_DEBUG
	    fprintf(stderr, "   (%d) * (log(%d/%d)) = %f\n", T->counter[a],
	       T->counter[a],
	       T->counter[3], s);
#endif*/
	    *likelihood += s;
	  }
      }

  /* descend on the tree on a depth-first manner */
  for (a = 0; a <= 2; a++)
    {
      if (T->sonsp[a] != NULL)
	calc_likelihood (T->sonsp[a], likelihood, k - 1);
    }
}


/*
 * bootstrap_memory: allocate space for likelihood, neighborhood and
 * log table data, and precalculate the log table. These actions are only
 * performed if the passed pointers are NULL.
 */
void
bootstrap_memory (double ***likelihood, pair_t ** neighborhoods, int max_nlen,
		  int ind_num)
{
  int i, j;
  double log3 = log (3);

  /* allocate space for likelihood table */
  if (*likelihood == NULL)
    {
      /* allocate space for table */
      *likelihood =
	(double **) malloc_or_fail ((max_nlen + 1) * sizeof (double *));
      (*likelihood)[0] =
	(double *) malloc_or_fail ((max_nlen + 1) * (max_nlen + 1) *
				   sizeof (double));
      /* adjust pointers */
      for (i = 1; i < max_nlen + 1; i++)
	(*likelihood)[i] = (*likelihood)[i - 1] + max_nlen + 1;
    }

  /* calculates logtable once */
  if (logtable == NULL)
    {
      /* allocate space for table */
      logtable =
	(double **) malloc_or_fail ((ind_num + 1) * sizeof (double *));
      logtable[0] =
	(double *) malloc_or_fail ((ind_num + 1) * (ind_num + 1) *
				   sizeof (double));
      /* adjust pointers */
      for (i = 1; i < ind_num + 1; i++)
	logtable[i] = logtable[i - 1] + ind_num;
      /* calculate log3(i/j) for every entry */
      for (i = 0; i < ind_num + 1; i++)
	for (j = 1; j < ind_num + 1; j++)
	  {
	    logtable[i][j] = log (((double) i) / ((double) j)) / log3;
	  }
    }

  /* allocates first space for neighborhood info */
  if (*neighborhoods == NULL)
    *neighborhoods = (pair_t *) malloc_or_fail (ind_num * sizeof (pair_t));
}
