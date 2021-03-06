/*
 Copyright (C) 2011 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id: points.cl 14305 2015-06-22 15:56:28Z dstrubbe $
*/

#include <cl_global.h>

__kernel void get_points(const int sp,
			 const int ep,
			 const int offset,
			 const int nst,
			 __global double const * restrict psi, const int ldpsi,
			 __global double * restrict points, const int ldpoints){
  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);
  const int sip = ip + sp - 1;

  if(ist < nst) points[ldpoints*ip + ist + offset] = psi[ldpsi*sip + ist];

}

__kernel void set_points(const int sp,
			 const int ep,
			 const int offset,
			 const int nst,
			 __global double const * restrict points, const int ldpoints,
			 __global double * restrict psi, const int ldpsi){
  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);
  const int sip = ip + sp - 1;

  if(ist < nst) psi[ldpsi*sip + ist] = points[ldpoints*ip + ist + offset];

}
/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
