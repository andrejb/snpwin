#
# common.mk - common definitions for all Makefiles
#
# Copyright (c) 2009 André Jucovsky Bianchi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

#-----------------------------------------------------------------------------
# This file contains common definitions for all Makefiles.
#
# See the README file for instructions.
#-----------------------------------------------------------------------------

RES_PREFIX := ../results
OBJ_PREFIX := ../common

#-----------------------------------------------------------------------------
# compiler definitions
#-----------------------------------------------------------------------------

CC := gcc
# for debugging:
CFLAGS := -O0 -lm -Wall -ansi -pedantic -g -std=c99
# for running:
#CFLAGS := -O3 -lm -Wall -ansi -pedantic

#-----------------------------------------------------------------------------
# main binaries
#-----------------------------------------------------------------------------

OBJ := $(subst .c,.o,$(wildcard $(OBJ_PREFIX)/*.c))

# rule for linking all binaries into one executable
$(BIN): $(BIN).o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

$(BIN).o: $(BIN).c
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJ_PREFIX)/%.o: %.c %.h

#-----------------------------------------------------------------------------
# other rules
#-----------------------------------------------------------------------------

all: target

clean:
	-\rm -f *.o

wipe: clean
	-\rm -f $(BIN)

edit:
	vim -p *.c *.h Makefile

indent:
	indent *.c
	if [ ! -z "`ls -1 | grep \\\.h`" ]; then indent *.h; fi
	rm *~

# targets that do not generate files
.PHONY: all clean wipe edit debug indent

