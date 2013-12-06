#
# Makefile - makefile for chi-square test
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
# This makefile has rules for (1) building and (2) running the chi-square test
# for testing the correlation between test and control windows.
#
# See the README file for instructions.
#-----------------------------------------------------------------------------

# (1) build rules

BIN := chi-square
include ../common/common.mk

# (2) rules for running the chi-square test

DATAFILE      = $(RES_PREFIX)/neighborhoods.txt
WINDOWSFILE   = $(RES_PREFIX)/windows-windows.txt
SOURCEFILE    = /var/tmp/data-ordered.txt
AFFECTIONFILE = ../data/affection.txt

run: do-chi-square

do-chi-square: chi-square
	./$< $(SOURCEFILE) $(WINDOWSFILE) $(AFFECTIONFILE) 

do-stats: stats
	cat $(WINDOWSFILE) | ./$<

