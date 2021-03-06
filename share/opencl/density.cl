/*
 Copyright (C) 2010 X. Andrade

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

 $Id: density.cl 14564 2015-09-13 19:26:50Z xavier $
*/

#include <cl_global.h>

__kernel void density_real(const int nst,
			   const int np,
			   const int offset,
			   __constant double * restrict weights,
			   const __global double * restrict psi, const int ldpsi,
			   __global double * restrict density){
  
  int ip  = get_global_id(0);
  if(ip >= np) return;

  double dd = 0.0;

  for(int ist = 0; ist < nst; ist++){
    double ff = psi[(ip<<ldpsi) + ist];
    dd += weights[ist]*ff*ff;
  }

  density[offset + ip] += dd;

}

__kernel void density_complex(const int nst,
			      const int np,
			      const int offset,
			      __constant double * restrict weights,
			      const __global double2 * restrict psi, const int ldpsi,
			      __global double * restrict density){
  
  int ip  = get_global_id(0);
  if(ip >= np) return;

  double dd = 0.0;

  for(int ist = 0; ist < nst; ist ++){
    double2 ff = psi[(ip<<ldpsi) + ist];
    ff = ff*ff;
    dd += weights[ist]*(ff.x + ff.y);
  }

  density[offset + ip] += dd;

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
