function H=HSMat(v,varflag)
%function H=HSMat(v,varflag)
%Input:  v = [v1;v2;d]  vi in R^(t(n)-n)  d in R^n  
%                  if varflag = complexCase
% v = [v1;d]  v1 in R^(t(n)-n)  d in R^n   
%                  if varflag = realCase 
% Output:  H  Hermitian if varflag == complexCase
%          H  Symmetric if varflag == realCase
       

if strcmp(varflag,'realCase')
    %n=round(sqrt(length(v)));  % FIX
    n = round( (-1+sqrt(1+8*length(v)) )/2); % inverse of triangular number function
    svind = find(triu(ones(n),1));  
    sqrt2=sqrt(2);

    tn = n*(n+1)/2; %number of variables in formulation
    tnm = tn-n;
    H=zeros(n);
    H(svind) = v(1:tnm)/sqrt2;
    H = (H+H');
    H(1:n+1:end) = v(tnm+1:end);
    H = (H+H')/2;
    
elseif strcmp(varflag,'complexCase')
    
    n=round(sqrt(length(v)));
    if n^2-length(v)>0
        disp('WARNING!!! catastrophic error; HMat index is wrong')
    end
    
    sqrt2=sqrt(2);
    svind = find(triu(ones(n),1)); 

    tn = n*(n+1)/2; %number of variables in formulation
    tnm = tn-n;
    H=zeros(n);
    H(svind) = (v(1:tnm) + 1i*v(tnm+1:2*tnm))/sqrt2;
    H = (H+H');
    H(1:n+1:end) = v(2*tnm+1:end);
    H = (H+H')/2;
end
end % of function
