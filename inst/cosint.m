## Copyright (C) 2006   Sissou   <sylvain.pelissier@gmail.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## COSINT compute the cosine integral function define by:
##
##                    Inf
##                   /
##       cosint(x) = | cos(t)/t dt
##                   /
##                   x
##
##See also : Ci, Si, sinint, expint, expint_Ei.

function y = cosint(z)
	if (nargin != 1)
		usage ("cosint(x)");
	endif
	y = Ci(z);
endfunction;