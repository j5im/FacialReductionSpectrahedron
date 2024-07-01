function c=CRvec(C,varflag)
%function c=CRvec(C)
% If varflag == realCase
%    INPUT:   C : n by n real matrix
%    OUTPUT:  c = C(:)
% If varflag == complesCase
%    INPUT:   C : n by n complex matrix
%    OUTPUT:  c = [real(C(:);imag(C(:)]
if strcmp(varflag,'realCase')
    c=C(:);
elseif strcmp(varflag,'complexCase')
    c=[real(C(:));imag(C(:))];
end
end  % end of function
