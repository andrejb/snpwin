#
# Makefile - general Makefile
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
# This makefile has rules for building all binaries.
#
# See the README file for instructions.
#-----------------------------------------------------------------------------

TARGETS = 1.neighborhoods \
          2.windows       \
          3.chi-square    \
          common

all:
	for i in $(TARGETS); do make -C $$i; done

clean:
	for i in $(TARGETS); do make wipe -C $$i; done
