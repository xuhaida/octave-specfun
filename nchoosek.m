## Copyright (C) 2001 Rolf Fabian, Paul Kienzle
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
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## c = nchoosek ( n,k )
##    Compute the binominal coefficient C(n,k) = n! / (k! (n-k)!)
##
## A = nchoosek ( v,k )
##    Generate all combinations of the elements of vector v
##    taken k at a time, one row per combination. The resulting
##    A has size [nchoosek(n,k), k], where n is the length of v.

##AUTHORS Rolf Fabian  <fabian@tu-cottbus.de>
##        Paul Kienzle <pkienzle@users.sf.net>

## XXX FIXME XXX shouldn't this be an alias for bincoeff?
## at the very least, compare implementations
function A = nchoosek(v,k)

  n = length(v);

  if n == 1
     A = round(exp( sum(log(k+1:v)) - sum(log(2:v-k)) ));

  elseif k == 0
    A = [];

  elseif k == 1
    A = v(:);

  elseif k == n
     A = v(:).';

  else

    m =  round(exp( sum(log(k:n-1)) - sum(log(2:n-k)) ));
    A = [ v(1)*ones( m,1 ), nchoosek(v(2:n), k-1) ; \
          nchoosek(v(2:n), k) ];
  endif

endfunction
