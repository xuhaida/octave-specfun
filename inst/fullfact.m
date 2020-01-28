## Copyright (C) 2012 - Michael Baudin
## Copyright (C) 2012 - Maria Christopoulou 
##
## This file must be used under the terms of the CeCILL.
## This source file is licensed as described in the file COPYING, which
## you should have received as part of this distribution.  The terms
## are also available at
## http:##www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function B = fullfact(levels)
    ## Full factorial design
    ##
    ## Calling Sequence
    ##   B = scidoe_fullfact(levels)
    ##
    ## Parameters
    ##   levels : a n-by-1 or 1-by-n matrix of doubles, integer value, positive, the number of levels for each factor j, for j=1,2,...,n
    ##   B : a m-by-n matrix of doubles, the experiments, where m=prod(levels). For j=1,2,...,n and i=1,2,...,m, and the entry B(i,j) is the level of the experiment #i for the variable #j. For a given variable j=1,2,...,n, the entries B(:,j) are in the set {1,2,...,level(j)}.
    ##
    ## Description
    ## Computes a full factorial design with prescribed number of 
    ## levels for each factor. 
    ## In other words, for each factor j, the number of levels is 
    ## levels(j), for j=1,2,...,n.
    ##
    ## Examples
    ## ## Create a full factorial design with :
    ## ## 2 levels for the first factor
    ## ## 3 levels for the second factor
    ## levels=[2 3]
    ## B=scidoe_fullfact(levels)
    ## ## Scale this design into [0,1]
    ## m=size(B,"r")
    ## C = (B-1)./(levels(ones(m,1),:)-1)
    ## ## Scale this design into [-1,1]
    ## D=2*C-1
    ## ## Plot this design
    ## scf();
    ## scidoe_plotcube(2);
    ## plot(D(:,1),D(:,2),"bo");
    ## xtitle("Full Factorial Design","X1","X2")
    ##
    ## ## Create a full factorial design with :
    ## ## 2 levels for the first factor
    ## ## 3 levels for the second factor
    ## ## 4 levels for the third factor
    ## levels = [2 3 4]
    ## B=scidoe_fullfact(levels)
    ## ## Scale this design into [0,1]
    ## m=size(B,"r")
    ## C = (B-1)./(levels(ones(m,1),:)-1)
    ## ## Scale this design into [-1,1]
    ## D=2*C-1
    ## ## Plot this design
    ## h = scf();
    ## param3d(D(:,1),D(:,2),D(:,3))
    ## h.children.children.mark_mode="on";
    ## h.children.children.line_mode="off";
    ## h.children.children.mark_size=1;
    ## scidoe_plotcube(3)
    ## xtitle("Full Factorial Design","X1","X2","X3")
    ##
    ## ## Print the number of experiments
    ## ## Use 3 levels for each parameter
    ## for n = 1 : 10
    ##   levels = 3*ones(n,1);
    ##   B=scidoe_fullfact(levels);
    ##   m = size(B,"r");
    ##   mprintf("n=%d, Num. Experiments=%d\n",..
    ##      n,m)
    ## end
    ##
    ## See also
    ## scidoe_plotcube
    ##
    ## Authors
    ## Copyright (C) 2012 - Michael Baudin
    ## Copyright (C) 2012 - Maria Christopoulou 

    ## Check number of input arguments
    ##[lhs,rhs] = argn();
    ##apifun_checkrhs("scidoe_fullfact",rhs,1);
    ##apifun_checklhs("scidoe_fullfact",lhs,1);
    ##
    ## Check type
    ## apifun_checktype("scidoe_fullfact",levels,"levels",1,["constant"]);
    ##
    ## Check size
    ## apifun_checkvector("scidoe_fullfact",levels,"levels",1);
    ##
    ## Check content
    ## apifun_checkflint("scidoe_fullfact",levels,"levels",1);
    ## apifun_checkgreq("scidoe_fullfact",levels,"levels",1,1);
    ##
    ## Proceed...
    args = {};
    for i = 1:numel(levels)
        args{end+1}=1:levels(i);
    end
    B = combvec ( args{1:end} )'
endfunction



%!test
%! levels = [2 3];
%! computed = fullfact(levels);
%! expected = [1 1;1 2;1 3;2 1;2 2;2 3];
%! assert(computed,expected);
