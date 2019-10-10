## Copyright (C) 2009 - Michael Baudin
## Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
##
## This file must be used under the terms of the CeCILL.
## This source file is licensed as described in the file COPYING, which
## you should have received as part of this distribution.  The terms
## are also available at
## http:##www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function c = combrepeat ( x , k )
    ##   Returns repeated combinations with replacement.
    ##
    ## Calling Sequence
    ##   c = specfun_combinerepeat ( x , k )
    ##
    ## Parameters
    ##   x : a m-by-n matrix, the matrix to produce combination from.
    ##   k : a 1-by-1 matrix of floating point integers, the number of repeted combinations, must be greater than 1.
    ##   c : a (k*m)-by-(n^k) matrix, same type as x
    ##
    ## Description
    ##   Uses a fast algorithm to produce repeated combinations of x with itself.
    ##   <itemizedlist>
    ##   <listitem>If k=1, then returns x.</listitem>
    ##   <listitem>If k=2, then returns combinations of x and x.</listitem>
    ##   <listitem>If k=3, then returns combinations of x and x and x.</listitem>
    ##   <listitem>etc...</listitem>
    ##   </itemizedlist>
    ##
    ##   For performance reasons, the combinations are stored column-by-column, that
    ##   is, the combination #k is in c(:,k), with k=1,2,...,n^k.
    ##
    ##   Can process x if x is double, strings boolean and integer (signed,unsigned,
    ##   8-bits, 16-bits, 32-bits).
    ##
    ##   The algorithm makes use of the Kronecker product for good performances.
    ##
    ## Examples
    ## ## Compute repeated combinations of x:
    ## x = [1 2 3];
    ## specfun_combinerepeat ( x , 1 )
    ## specfun_combinerepeat ( x , 2 )
    ## specfun_combinerepeat ( x , 3 )
    ## specfun_combinerepeat ( x , 4 )
    ##
    ## ## Compare to specfun_combine
    ## ## Same as k=2
    ## specfun_combine ( x , x )
    ## ## Same as k=3
    ## specfun_combine ( x , x , x )
    ## ## Same as k=4
    ## specfun_combine ( x , x , x , x )
    ##
    ## ## Repeated combinations of booleans
    ## computed = specfun_combinerepeat ( [%t %f] , 2 )
    ## ## Repeated combinations of strings
    ## computed = specfun_combinerepeat ( ["A" "C" "T" "G"] , 2 )
    ## ## Repeated combinations of integers
    ## computed = specfun_combinerepeat ( uint8(1:3) , 2 )
    ##
    ## ## Compare combinerepeat, perms and subset
    ## ## Scilab/perms compute permutations without replacement
    ## perms ( 1:3 )
    ## ## specfun_combinerepeat compute combinations with replacement
    ## specfun_combinerepeat(1:3,3)'
    ## ## specfun_subset compute subsets with k elements
    ## specfun_subset(1:3,3)
    ##
    ## Authors
    ## Copyright (C) 2009 - Michael Baudin
    ## Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
    ##

    ##[lhs,rhs]=argn()
    ## apifun_checkrhs ( "specfun_combinerepeat" , rhs , 2 )
  ##  apifun_checklhs ( "specfun_combinerepeat" , lhs , 0:1 )
    ##
    ## Check Type
   ## apifun_checktype ( "specfun_combinerepeat" , x , "x" , 1 , ["constant" "string" "boolean" "int8" "uint8"  "int16" "uint16" "int32" "uint32"] )
   ## apifun_checktype ( "specfun_combinerepeat" , k , "k" , 2 , "constant" )
    ##
    ## Check size
   ## apifun_checkscalar ( "specfun_combinerepeat" , k , "k" , 2 )
    ##
    ## Check content
   ## apifun_checkgreq ( "specfun_combinerepeat" , k , "k" , 2 , 1 )
   ## apifun_checkflint ( "specfun_combinerepeat" , k , "k" , 2 )
    ## TODO : use apifun_checkflint when ready
    ##
    ## Proceed...
  c = x;
    if ( k == 1 ) 
      return
    end
    for j = 2 : k
      c = combvec( c , x );
    end
endfunction

%!test

## Test 2 combinations with one row vector
%! x = [1 2 3];
%! computed = combrepeat ( x , 2 );
%! expected = [
%! 1,1,1,2,2,2,3,3,3;
%! 1,2,3,1,2,3,1,2,3
%! ];
%! assert( computed , expected );

##Test 1 combinations with one row vector
%! x = [1 2 3];
%! computed = combrepeat ( x , 1 );
%! expected = x;
%! assert( computed , expected );

## Test 3 combinations with one row vector
%! x = [1 2 3];
%! computed = combrepeat ( x , 3 );
%! expected = [
%!   1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3;
%!   1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3;
%!   1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3
%! ];
%! assert( computed , expected );

## Repeated combinations of booleans
%! computed = combrepeat ( [true  false] , 2 );
%! expected = [true, true , false,false;true , false, true , false];
%! assert( computed , expected );
##
## Repeated combinations of strings
%! computed = combrepeat ( ["A" "C" "T" "G"] , 2 );
%! expected = [
%! "A","A","A","A","C","C","C","C","T","T","T","T","G","G","G","G";
%! "A","C","T","G","A","C","T","G","A","C","T","G","A","C","T","G"
%! ];
%! assert( computed , expected );
##
## Repeated combinations of integers
##%! computed = combrepeat ( uint8(1:3) , 2 );
##%! expected = uint8([
##%! 1,1,1,2,2,2,3,3,3;
##%! 1,2,3,1,2,3,1,2,3
##%! ]);
##%! assert( computed , expected );


