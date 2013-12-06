/*
 * tree.h - definitions for tree.c
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

#ifndef _tree_h
#define _tree_h

/*
 * tree_node: holds (1) counters for the number of times each symbol appeared
 * in a specific position of a string, and (2) pointers for next possible
 * symbols.
 */
typedef struct tree_struct
{
  int counter[4];
  struct tree_struct *sonsp[3];
  struct tree_struct *next;
} tree_node;

/* pointer for a tree node */
typedef tree_node *tree_type;

/* public functions */
tree_type make_empty_tree ();
void print_tree (tree_type T, int level, int label);
void insert_string (tree_type T, char **snpxind, int ind, int snp,
		    int snps_left, int a, int jump);
void free_tree (tree_type T);
int snptoi (char snp);

#endif /*_tree_h */
