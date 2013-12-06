/*
 * read_data.h - definitions for read_data.c
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

#ifndef _read_data_h
#define _read_data_h

#define STRLEN_MAX 20000	/* technical, buffer upper limit */
			 /* mustn't exceed <limits.h> INT_MAX */
			 /* in v1.2 and v2.0 longest string = 5217 */

/* actualy we also assume sum of all string length to fit into an int
   therefore if sizeof(int)=2 expect to fail on long data-set files */


void read_data (char ***snpxind,
		char **data, FILE * stream, int snp_num, int ind_num);

FILE *open_file (char *fname, char *mode);

#define SNP_NUM    500180
#define IND_NUM    2062

#endif /* _read_data_h */
