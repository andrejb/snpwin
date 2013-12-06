/*
 * tree.c - tree primitives for counting strings
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

#include <math.h>
/*#include <unistd.h>*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tree.h"

/* local */

/* global node pool */
tree_node *pool = NULL;

/*
 * reset_node: reset node's counters and pointers.
 */
void
reset_node (tree_node * nodep)
{
  nodep->counter[0] = 0;
  nodep->counter[1] = 0;
  nodep->counter[2] = 0;
  nodep->counter[3] = 0;
  nodep->sonsp[0] = NULL;
  nodep->sonsp[1] = NULL;
  nodep->sonsp[2] = NULL;
}

/*
 * tree_node: returns a new node for inserting SNP values on the tree.
 */
tree_node *
make_tree_node ()
{
  static int pool_size = 0;
  tree_node *new_pool;
  tree_node *nodep;
  int i, s;
  /* we ran out of nodes! */
  if (pool == NULL)
    {
      s = (pool_size == 0 ? 2 : pool_size);
      new_pool = (tree_node *) malloc (s * sizeof (tree_node));
      if (new_pool == NULL)
	{
	  fprintf (stderr, "failed to allocate tree node\n");
	  exit (1);
	}
      for (i = 0; i < s; i++)
	{
	  reset_node (&(new_pool[i]));
	  new_pool[i].next = ((i == s - 1) ? NULL : &(new_pool[i + 1]));
	}
      pool = new_pool;
      pool_size += s;
    }
  nodep = pool;
  pool = pool->next;
  nodep->next = NULL;
  return nodep;
}

/*
 * free_tree: return every node on a tree to the node pool.
 */
void
free_tree (tree_type Tp)
{
  int i;
  if (Tp != NULL)
    {
      Tp->next = pool;
      pool = Tp;
      for (i = 0; i <= 2; i++)
	free_tree (Tp->sonsp[i]);
      reset_node (Tp);
    }
  return;
}


/* public */

/*
 * snptoi: convert char values read from SNP input data to integers.
 */
int
snptoi (char snp)
{
  if (snp == '3')
    return -1;
  return (int) (snp - 48);
}

/*
 * insert_string: inserts a string of SNP values into the tree, possibly
 * skipping one of the values (this is needed for the calculation of the
 * likelihood of a neighborhood, where we do not consider the central SNP).
 */
void
insert_string (tree_type T,
                char **snpxind,
                int ind,
                int snp,
                int snps_left,
                int a,
                int skip)
{
  int son;

  /* verify is SNP value is valid */
  if (a != 0 && a != 1 && a != 2)
    {
      fprintf (stderr, "Error: symbol not recognized: %d.\n", a);
      exit (1);
    }

  /* count the SNP value as occurring once more */
  T->counter[a] += 1;
  T->counter[3] += 1;
  /*printf("%d SNPs missing and will skip the %d-th.\n", snps_left, skip); */

  /* insert the rest of symbols of the string (if any) */
  if (snps_left > 0)
    {
      /* skip one SNP if this is the case */
      if (skip == 0)
	snp++;

      son = snptoi (snpxind[snp + 1][ind]);
      /*printf("inserindo o filho %d\n", son); */
      if (T->sonsp[son] == NULL)
	T->sonsp[son] = make_tree_node ();

      /* recursivelly insert the rest of the string */
      insert_string (T->sonsp[son], snpxind, ind, snp + 1, snps_left - 1, a,
		     skip - 1);
    }
}

/*
 * make_empty_tree: returns a fresh new tree.
 */
tree_type
make_empty_tree ()
{
  return (make_tree_node ());
}

/*
 * print_tree: print the whole tree (mainly for debugging).
 */
void
print_tree (tree_type T, int level, int label)
{
  int i;

  if (level == 0)
    printf ("root\n");

  for (i = 0; i < level; i++)
    printf ("    ");
  printf ("[%d, %d (%d, %d, %d)] -> [ ", label, T->counter[3], T->counter[0],
	  T->counter[1], T->counter[2]);
  for (i = 0; i <= 2; i++)
    if ((T->sonsp)[i] != NULL)
      printf ("%d ", i);
  printf ("]\n");
  ++level;
  for (i = 0; i <= 2; i++)
    if ((T->sonsp)[i] != NULL)
      print_tree ((T->sonsp)[i], level, i);
}
