## Copyright (C) 2001 Paul Kienzle
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

## A = perms(v)
##    generate all permutations of v, one row per permutation.
##    The resulting A has size n! x n, where n is the length of v
##    so keep v small!
function A = perms(v)
  n = length(v);
  if (n == 1)
    A = v;
  else
    B = perms(v(1:n-1));
    Bidx = 1:size(B,1);
    A = v(n) * ones(prod(2:n), n);
    A (Bidx, 1:n-1) = B;
    k = size(B,1);
    for i = n-1:-1:2
      A (k+Bidx, 1:i-1) = B(Bidx, 1:i-1);
      A (k+Bidx, i+1:n) = B(Bidx, i:n-1);
      k = k + size(B,1);
    end
    A (k+Bidx, 2:n) = B;
  endif
endfunction
