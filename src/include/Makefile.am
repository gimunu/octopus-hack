## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: Makefile.am 15058 2016-01-14 16:50:31Z xavier $

noinst_HEADERS =                  \
	global.h                  \
	integer.F90               \
	real.F90                  \
	complex.F90               \
	real_single.F90           \
	complex.F90               \
	complex_single.F90        \
	undef.F90                 \
	vectors.h                 \
	fortran_types.h

include $(srcdir)/../common-rules.make

global.h: config_F90.h options.h defaults.h

config_F90.h: $(top_builddir)/config.h
	grep '^#define' $^ > config_F90.h

options.h: RUN_SCRIPT

defaults.h: RUN_SCRIPT

RUN_SCRIPT:
	$(top_srcdir)/build/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

clean-local:
	rm -f config_F90.h options.h defaults.h

.PHONY: RUN_SCRIPT
