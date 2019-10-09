
## Copyright (C) 2009 - Michael Baudin
## Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
## Copyright (C) 2019 - Haida Xu

## This file must be used under the terms of the CeCILL.
## This source file is licensed as described in the file COPYING, which
## you should have received as part of this distribution.  The terms
## are also available at
## http:##www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function c = combvec ( varargin )
    ##   Returns the all the combinations of the given vectors.
    ##
    ## Calling Sequence
    ##   c = combvec ( a1 , a2 )
    ##   c = combvec ( a1 , a2 , a3 )
    ##   c = combvec ( a1 , a2 , a3 , ... )
    ##
    ## Parameters
    ##   a1 : a m1-by-n1 matrix
    ##   a2 : a m2-by-n2 matrix
    ##   a3 : a m3-by-n3 matrix
    ##   c : a m-by-n matrix, with m=m1+m2+m3+... and n=n1*n2*n3*...
    ##
    ## Description
    ##   Uses a fast algorithm to produce combinations.
    ##
    ##   For performance reasons, the combinations are stored column-by-column, that
    ##   is, the combination #k is in c(:,k), with k=1,2,...,n1*n2*n3*....
    ##   The rows 1 to m1 contains elements from a1.
    ##   The rows m1+1 to m1+m2 contains elements from a2.
    ##   etc...
    ##
    ##   Can process matrices of doubles, strings, booleans and all integers (signed,unsigned,
    ##   8-bits, 16-bits, 32-bits,64-bits).
    ##
    ##   The algorithm makes use of the Kronecker product for good performances.
    ##
    ##   It works, for example, if:
    ##   <itemizedlist>
    ##   <listitem>all arguments are row vectors, with same number of entries,</listitem>
    ##   <listitem>all arguments are column vectors, with same number of entries,</listitem>
    ##   <listitem>all arguments are matrices, with same number of rows and columns,</listitem>
    ##   <listitem>a row vector and a scalar,</listitem>
    ##   <listitem>etc...</listitem>
    ##   </itemizedlist>
    ##
    ## Examples
    ## ## Compute all combinations of x and y: vectors
    ## x = [1 2 3];
    ## y = [4 5 6];
    ## specfun_combine ( x , y )
    ## expected = [
    ##   1 1 1 2 2 2 3 3 3
    ##   4 5 6 4 5 6 4 5 6
    ## ];
    ##
    ## ## Compute all combinations of x and y and z: vectors
    ## x = [1 2 3];
    ## y = [4 5 6];
    ## z = [7 8 9];
    ## specfun_combine ( x , y , z )
    ## expected = [
    ##   1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3
    ##   4 4 4 5 5 5 6 6 6 4 4 4 5 5 5 6 6 6 4 4 4 5 5 5 6 6 6
    ##   7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9
    ## ];
    ##
    ## ## Compute all combinations of x and y: matrices
    ## x = [1 2;3 4];
    ## y = [5 6;7 8];
    ## specfun_combine ( x , y )
    ## expected = [
    ##   1 1 2 2
    ##   3 3 4 4
    ##   5 6 5 6
    ##   7 8 7 8
    ## ];
    ##
    ## ## Combine random matrices of random sizes.
    ## ## Shows that matrices of any dimensions can be combined.
    ## m = grand(1,2,"uin",1,5);
    ## n = grand(1,2,"uin",1,5);
    ## x = grand(m(1),n(1),"uin",1,m(1)*n(1));
    ## y = grand(m(2),n(2),"uin",1,m(2)*n(2));
    ## c = specfun_combine ( x , y );
    ## and(size(c) == [m(1)+m(2) n(1)*n(2)])
    ##
    ## ## Indirectly produce combinations of characters
    ## k = specfun_combine ( 1:2 , 1:2 )
    ## m1=["a" "b"];
    ## m2=["c" "d"];
    ## c = [m1(k(1,:));m2(k(2,:))]
    ##
    ## ## Directly combine strings
    ## x = ["a" "b" "c"];
    ## y = ["d" "e"];
    ## z = ["f" "g" "h"];
    ## computed = specfun_combine ( x , y , z )
    ##
    ## ## Produces combinations of booleans
    ## c = specfun_combine ( [%t %f] , [%t %f] , [%t %f] )
    ##
    ## ## Combine 2 DNA genes
    ## c = specfun_combine ( ["A" "C" "G" "T"]  , ["A" "C" "G" "T"] )
    ##
    ## ## Produces combinations of integers
    ## c = specfun_combine(uint8(1:4),uint8(1:3))
    ##
    ## Authors
    ## Copyright (C) 2009 - Michael Baudin
    ## Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
    ## Copyright (C) 2019 - Haida Xu

    ##TODO: to be delete
##    function c = combine2argsDouble_old ( X , Y )
        ## Returns all combinations of the two row vectors X and Y
##        cx=size(X,"c")
##        cy=size(Y,"c")
##        c=[X .*. ones(1,cy); ones(1,cx) .*. Y]
##    endfunction


##    narginchk( 1,32 )
##    nargoutchk( 0,1 )

  rhs=nargin();
##    lhs=nargout()
    
    ##
  ## Check that all arguments are constants or strings
  for k=1:rhs
    
    validateattributes( varargin{k}, {"double" "char" "boolean" "int8" "uint8"  "int16" "uint16" "int32" "uint32"},{'nonempty'} );
  end
    ##
    ## Check that all arguments have the same type
  type1=class(varargin{1});
  for k=2:rhs
    validateattributes(varargin{k},{type1},{'nonempty'});
  end
  
    ##
    if ( rhs == 1 ) then
      c = varargin(1);
        return
    end
    if ( strcmp(type1,"double") ) 
      ##
  ## Combine doubles
  
      c = combine2argsDouble ( varargin{1} , varargin{2} );
      for i = 3:rhs
	c = combine2argsDouble ( c , varargin{i} );
      end
    else
      ##
      ## Combine strings
      c = combine2argsOther ( varargin{1} , varargin{2} );
      for i = 3:rhs
        c = combine2argsOther ( c , varargin{i} );
      end
    end
endfunction

function c = combine2argsDouble ( x , y )
  
   ## Returns all combinations of the two row vectors x and y
  cx=size(x,2);
  cy=size(y,2);
   ##
   ##
  i1 = 1:cx;
  j1 = i1(ones(cy,1),:);
  j1 = reshape(j1,1,cx*cy);
   ##
  i2 = 1:cy;
  j2 = i2(ones(cx,1),:);
  j2 = reshape(j2',1,cx*cy);
   ##
   c = [
      x(:,j1)
      y(:,j2)
   ];
endfunction

function c = combine2argsOther ( X , Y )
      ##
      ## Compute the size of the args
  m(1) = size(X,1);
  m(2) = size(Y,1);
  n(1) = size(X,2);
  n(2) = size(Y,2);
      ##
      ## Combine the two first arguments
      ##
      ## Create and combine a matrix of indices
  matrix1 = reshape(1:m(1)*n(1),m(1),n(1));
  matrix2 = reshape(1:m(2)*n(2),m(2),n(2));
  cm = combine2argsDouble ( matrix1 , matrix2 );
  
      ##
      ## Get the actual strings to combine
  j = (1:n(1)*n(2));
      for i = 1:m(1)
        d = X(cm(i,j));
        c(i,j) = d(:)';
      end
      for i = 1:m(2)
        d = Y(cm(i+m(1),j))';
        c(i+m(1),j) = d(:)';
      end
 endfunction

 


%!test
%! x = [1 2 3];
%! y = [4 5 6];
%! expected = [
%!   1 1 1 2 2 2 3 3 3
%!   4 5 6 4 5 6 4 5 6
%!  ];
%! assert(combvec(x,y),expected);

## Compute all combinations of x and y and z: vectors


%! z = [7 8 9];
%! expected = [
%!   1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3
%!   4 4 4 5 5 5 6 6 6 4 4 4 5 5 5 6 6 6 4 4 4 5 5 5 6 6 6
%!   7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9 7 8 9
%! ];
%! assert(combvec(x,y,z),expected);

## ## Indirectly produce combinations of characters
## m1=["a" "b"];
## m2=["c" "d"];
## expacted=[

## ## Compute all combinations of x and y: matrices
%! x = [1 2;3 4];
%! y = [5 6;7 8];

 
%! expected = [
%!   1 1 2 2
%!   3 3 4 4
%!   5 6 5 6
%!   7 8 7 8
%! ];
%! assert(combvec(x,y),expected);

    ## ## Combine 2 DNA genes
%! x= ["A" "C" "G" "T"];
%! y= ["A" "C" "G" "T"];
%! expected=[
%! "AAAACCCCGGGGTTTT"
%! "ACGTACGTACGTACGT"
%! ];
%! assert(combvec(x,y),expected);


## ## Directly combine strings
%! x = ["a" "b" "c"];
%! y = ["d" "e"];
%! z = ["f" "g" "h"];
%! expected=[
%!   "aaaaaabbbbbbcccccc"
%!   "dddeeedddeeedddeee"
%!   "fghfghfghfghfghfgh"
%! ];
%! assert(combvec(x,y,z),expected);


    
