## Copyright (C) 2000  Kai Habel
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

## -*- texinfo -*-
## @deftypefn {Function File} {@var{L} =}legendre (@var{n},@var{X})
## Legendre Function of degree n and order m
## where all values for m = 0..@var{n} are returned.
## @var{n} must be a scalar in the range [0..255].
## The return value has one dimension more than @var{x}.
##
## @example
## The Legendre Function of degree n and order m
##
## @group
##  m        m       2  m/2   d^m
## P(x) = (-1) * (1-x  )    * ----  P (x)
##  n                         dx^m   n
## @end group
##
## with:
## Legendre Polynom of degree n
##
## @group
##           1     d^n   2    n
## P (x) = ------ [----(x - 1)  ] 
##  n      2^n n!  dx^n
## @end group
##
## legendre(3,[-1.0 -0.9 -0.8]) returns the matrix
##
## @group
##  x  |   -1.0   |   -0.9   |  -0.8
## ------------------------------------
## m=0 | -1.00000 | -0.47250 | -0.08000
## m=1 |  0.00000 | -1.99420 | -1.98000
## m=2 |  0.00000 | -2.56500 | -4.32000
## m=3 |  0.00000 | -1.24229 | -3.24000 
## @end group
## @end example
## @end deftypefn

## Author:	Kai Habel <kai.habel@gmx.de>

function L = legendre (n,x)
  try dfi = do_fortran_indexing;
  catch dfi = 0;
  end
  try wfi = warn_fortran_indexing;
  catch wfi = 0;
  end
  unwind_protect
    do_fortran_indexing = 1;
    warn_fortran_indexing = 0;
    warning ("legendre is unstable for higher orders");
    if (nargin != 2)
      usage ("legendre(n,x)");
    endif

    if !(is_scalar (n) && n >= 0 && n < 256 && n == fix (n))
      error ("n must be a scalar n=[0...255]");
    endif

    if !(is_vector (x) && all (x >= -1 && x <= 1))
      error ("x must be vector in range -1<=x<=1");
    endif
    if (n == 0)
      L = ones (size (x));
    elseif (n == 1)
      L = [x; -x .* sqrt(1 - x .^ 2)];
    else
      i = 0:n;
      a = (-1) .^ i .* bincoeff (n, i);
      p = [a; zeros(size (a))];
      p = p(:);
      p(length (p)) = [];
      #p contains the polynom (x^2-1)^n

      #now create a vector with 1/(2.^n*n!)*(d/dx).^n
      d = [((n + rem(n, 2)):-1:(rem (n, 2) + 1)); 2 * ones(fix (n / 2), n)];
      d = cumsum (d);
      d = [fliplr(prod (d'))];
      d = [d; zeros(1, length (d))];
      d = d(1:n + 1) ./ (2 ^ n *prod (1:n));

      Lp = d' .* p(1:length (d));
      #Lp contains the Legendre Polynom of degree n

      # now create a polynom matrix with d/dx^m for m=0..n
      d2 = flipud (triu (ones (n)));
      d2 = cumsum (d2);
      d2 = fliplr (cumprod (flipud (d2)));
      d3 = fliplr (triu (ones (n + 1)));
      d3(2:n + 1, 1:n) = d2;

      # multiply for each m (d/dx)^m with Lp(n,x)
      # and evaluate at x
      Y = zeros(n + 1, length (x));
      [dr, dc] = size (d3);
      for m = 0:dr - 1
        Y(m + 1, :) = polyval (d3(m + 1, 1:(dc - m)) .* Lp(1:(dc - m))', x);
      endfor

      # calculate (-1)^m*(1-x^2)^(m/2)	for m=0..n at x
      # and multiply with (d/dx)^m(Pnx)
      l = length (x);
      X = kron ((1 - x(:) .^ 2)', ones (n + 1, 1));
      M = kron ((0:n)', ones (1, l));
      L = X .^ (M / 2) .* (-1) .^ M .* Y;
    endif
  unwind_protect_cleanup
    do_fortran_indexing = dfi;
    warn_fortran_indexing = wfi;
  end_unwind_protect
endfunction
