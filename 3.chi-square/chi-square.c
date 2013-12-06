/*
 * arquivo.c - descricao do arquivo.c
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
#include <unistd.h>
#include <getopt.h>

#include "../common/tree.h"
#include "../common/util.h"
#include "../common/read_data.h"

#define OUTPUT_FILE "../results/chi-square-chi-square.txt"

#define GROUP_DEFAULT_VALUE 5
#define OPT_STRING "+g::G"

#define max(x,y) (x > y ? x : y)

/* function prototypes */
int read_windows (char *filename, int *windows);
int read_affection (char *filename, int *affection);
double chi_square (int start, int wlen, char **snpxind,
		   pair_t * neighborhoods, int *affection, FILE * fp, int
                   group);
void calc_chi_square (tree_type T, int count_cases, int count_control, int lfcounter[2], int k,
		      double *chi, int *df, int group);
void usage(int exit_status, char *program_name);

/* main: entry point for SNP windows chi-square statistic calculation.
 *
 * TODO
 */
int
main (int argc, char *argv[])
{

  FILE *infile, *outfile;	/* handlers for input and output files */
  char **snpxind;
  char *data;
  int windows[SNP_NUM], num_windows, snp, i;
  int affection[IND_NUM];
  pair_t neighborhoods[IND_NUM];
  int opt, group = 0;
  /*int superposition[SNP_NUM]; */

  /* check number of command-line arguments */
  if (argc < 4 || argc > 5)
    {
      printf ("%s: invalid number of args (%d).\n", argv[0], argc-1);
      usage(1, argv[0]);
    }

  /* get grouping option */
  if ((opt = getopt(argc, argv, OPT_STRING)) != -1) {
    if (opt == 'g')
      if (optarg != 0)
        group = atoi(optarg);
      else
        group = GROUP_DEFAULT_VALUE;
    else if (opt == 'G')
      group = 0;
    else {
      usage(1, argv[0]);
    }
    if (getopt(argc, argv, OPT_STRING) != -1) {
      printf ("%s: can't have -g and -G simultaneously.\n", argv[0]);
      usage(1, argv[0]);
    }
    fprintf(stderr, "Will group strings with expected frequency less than %d.\n", group);
  }

  /* open file for output */
  outfile = open_file (OUTPUT_FILE, "w");

  /* read input data */
  infile = open_file (argv[optind], "r");
  read_data (&snpxind, &data, infile, SNP_NUM, IND_NUM);
  fclose (infile);

  /* read previously calculated data */
  num_windows = read_windows (argv[optind+1], windows);
  read_affection (argv[optind+2], affection);

  /* calculate chi-square for each window */
  snp = 0;
  for (i = 0; i < num_windows; i++)
    {
      /* limit neighborhoods based on NA values */
      limit_neighborhoods (neighborhoods, snpxind, snp, IND_NUM, 0,
			   windows[i] - 1);
      /* perform chi-square test */
      chi_square (snp, windows[i], snpxind, neighborhoods, affection,
		  outfile, group);
      /* advance to next window */
      snp += windows[i];
    }

  return 0;
}



/*
 * read_windows:
 */
int
read_windows (char *filename, int *windows)
{
  FILE *infile;
  char buf[STRLEN_MAX];
  int start, end, count;

  /* open the file with window info */
  infile = open_file (filename, "r");

  /* read info from file */
  count = 0;
  while ((fgets (buf, STRLEN_MAX, infile) != NULL) && buf[0] != '\0')
    {
      sscanf (buf, "%d %d", &start, &end);
      windows[count] = end - start + 1;
      count++;
    }

  return count;
}

/*
 * read_affection:
 */
int
read_affection (char *filename, int *affection)
{
  FILE *infile;
  char buf[STRLEN_MAX];
  int indicator, count;

  /* open the file with affection info */
  infile = open_file (filename, "r");

  /* read info from file */
  count = 0;
  while ((fgets (buf, STRLEN_MAX, infile) != NULL) && buf[0] != '\0')
    {
      sscanf (buf, "%d", &indicator);
      affection[count] = indicator;
      count++;
    }

  return count;
}

/*
 * chi_square:
 */
double
chi_square (int start,		/* beginning of the window */
	    int wlen,		/* window length */
	    char **snpxind,	/* SNP data */
	    pair_t * neighborhoods,	/* neighborhoods data */
	    int *affection,	/* affection data */
	    FILE * fp,		/* output file pointer */
            int group) /* group outcomes with less than this amount of expected frequency */
{
  int i, count_cases, count_control, df, lfcounter[2];
  tree_type T;
  double chi;
  double EK, EL, TK, TL;

  /* initialize variables */
  T = make_empty_tree ();
  count_cases = 0;
  count_control = 0;
  lfcounter[0] = 0;
  lfcounter[1] = 0;

  /* iterate for every individual in the sample */
  for (i = 0; i < IND_NUM; i++)
    {
      /* test if this individual's window is valid: */
      if (snptoi (snpxind[start][i]) != -1	/* first symbol is not NA */
	  && neighborhoods[i].r == wlen - 1)	/* rest of the symbols are not NA */
	{
	  /* insert window and count individual as case or control */
	  insert_string (T, snpxind, i, start - 1, wlen, affection[i], -1);
	  if (affection[i] == 1)
	    count_cases++;
	  else
	    count_control++;
	}
    }

  /* calculate chi-square statistic for the tree */
  chi = 0;
  df = 0;
  calc_chi_square (T, count_cases, count_control, lfcounter, wlen, &chi, &df, group);
  free_tree (T);

  if (group > 0) {
    /* calculate chi-square for least-frequent words */
    EK =
      ((double) ((lfcounter[0] + lfcounter[1]) * count_cases)) /
      (count_cases + count_control);
    EL =
      ((double) ((T->counter[0] + T->counter[1]) * count_control)) /
      (count_cases + count_control);
    if (EL != 0)
      TL = ((double) (lfcounter[0] - EL) * (lfcounter[0] - EL)) / EL;
    else
      TL = 0;
    if (EK != 0)
      TK = ((double) (lfcounter[1] - EK) * (lfcounter[1] - EK)) / EK;
    else
      TK = 0;

    /* sum chi-square of least-frequent words to the overall statistic */
    if (TK > 0 || TL > 0) {
      chi += TK + TL;
      df += 1;
    }
  }
  
  /* print the results */
  fprintf (fp, "%f %d\n", chi, (int) max(df - 1, 1));

  return chi;
}


/*
 * calc_chi_square:
 */
void
calc_chi_square (tree_type T,
		 int count_cases,
		 int count_control, int lfcounter[2], int k, double *chi, int
                 *df, int group)
{
  int a;
  double EK,			/* expected number of case individuals */
    EL,				/* expected number of control individuals */
    TK,				/* total number of case individuals */
    TL;				/* total number of control individuals */

  /* if we are on a leaf, include the word of the branch on the statistics */
  if (k == 0)
    {
      /* calculate expected number of case individuals */
      EK =
	((double) ((T->counter[0] + T->counter[1]) * count_cases)) /
	(count_cases + count_control);
      /* calculate expected number of control individuals */
      EL =
	((double) ((T->counter[0] + T->counter[1]) * count_control)) /
	(count_cases + count_control);

      /* verifies if this string has enough frequency... */
      if (fmin(EK, EL) >= group) {
        /* calculate total number of control individuals */
        if (EL != 0)
          TL = ((double) (T->counter[0] - EL) * (T->counter[0] - EL)) / EL;
        else
          TL = 0;
        /* calculate total number of case individuals */
        if (EK != 0)
          TK = ((double) (T->counter[1] - EK) * (T->counter[1] - EK)) / EK;
        else
          TK = 0;

        /* sum parcel to chi-square statistics */
        *chi += TK + TL;
        *df += 1;
      }
      /* ... otherwise, group this word with others with that also have low frequency. */
      else {
        lfcounter[0] += T->counter[0];
        lfcounter[1] += T->counter[1];
      }
    }

  /* if not on a leaf yet, go down on the branch */
  else
    for (a = 0; a <= 2; a++)
      if (T->sonsp[a] != NULL)
	calc_chi_square (T->sonsp[a], count_cases, count_control, lfcounter, k - 1, chi,
			 df, group);
}

/* usage: print info about how to use the program and exit with given status.
 */
void usage(int exit_status, char *program)
{
  printf("Usage: %s [-g [n] | -G] <source_file> <windows_file> <affection_file>\n", program);
  printf("  -G    do not group different words (default).\n");
  printf("  -g    group words with expected frequency less than n (5 if no frequency is\n");
  printf("        informed.\n");
  exit(exit_status);
}
