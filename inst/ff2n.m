## Copyright (C) 2012 - Michael Baudin
## Copyright (C) 2012 - Maria Christopoulou 
##
## This file must be used under the terms of the CeCILL.
## This source file is licensed as described in the file COPYING, which
## you should have received as part of this distribution.  The terms
## are also available at
## http:##www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function B = ff2n(n)
    ## Full factorial design with 2 levels
    ##
    ## Calling Sequence
    ##   B=scidoe_ff2n(n)
    ##
    ## Parameters
    ##   n : a 1-by-1 matrix of doubles, integer value, positive, the number of levels for each factor
    ##   B : a m-by-n matrix of doubles, the experiments in the range [0,1], where m=n*n. 
    ##
    ## Description
    ## Computes a full factorial design with 2 levels for each factor.
    ##
    ## Examples
	## ## Create a full factorial design with :
	## ## 2 levels for the first factor
	## ## 2 levels for the second factor
    ## B=scidoe_ff2n(2)
	## ## Scale into [-1,1]
	## C=2*B-1;
	## ## Plot this design
	## scf();
	## scidoe_plotcube(2);
	## plot(C(:,1),C(:,2),"bo");
	## xtitle("Full Factorial Design","X1","X2")
	##
	## ## For three factors
    ## B=scidoe_ff2n(3)
	## ## Scale into [-1,1]
	## C=2*B-1;
	## ## Plot this design
	## h = scf();
	## param3d(C(:,1),C(:,2),C(:,3))
	## h.children.children.mark_mode="on";
	## h.children.children.line_mode="off";
	## h.children.children.mark_size=1;
	## scidoe_plotcube(3)
	## xtitle("Full Factorial Design","X1","X2","X3")
	##
	## ## For four factors
    ## B=scidoe_ff2n(4)
	##
	## ## Print the number of experiments
	## for n = 1 : 10
	##   B=scidoe_ff2n(n);
	##   m = size(B,"r");
	##   mprintf("n=%d, Num. Experiments=%d\n",..
	##      n,m)
	## end
    ##
    ## Authors
	## Copyright (C) 2012 - Michael Baudin
	## Copyright (C) 2012 - Maria Christopoulou 

	##
    ## Check number of input arguments
    ## [lhs,rhs] = argn();
    ## apifun_checkrhs("scidoe_ff2n",rhs,1);
    ## apifun_checklhs("scidoe_ff2n",lhs,1);
	##
    ## Check type
    ## apifun_checktype("scidoe_ff2n",n,"n",1,["constant"]);
	##
    ## Check size
    ## apifun_checkscalar("scidoe_ff2n",n,"n",1);
	##
    ## Check content
    ## apifun_checkflint("scidoe_ff2n",n,"n",1);
    ## apifun_checkgreq("scidoe_ff2n",n,"n",1,1);
	##
	## Proceed...
    B = combrepeat ( [0,1] , n)'
endfunction

%!test
%! n=3;
%! computed = ff2n(n);
%! expected = [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];
%! assert(computed,expected);
