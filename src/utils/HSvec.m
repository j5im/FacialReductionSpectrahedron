function v=HSvec(H,varflag)
%function v=HSvec(H,varflag)
%INPUT:   H n by n Hermitian if varflag == complexCase
%         H n by n symmetric if varflag == realCase
%OUTPUT:  v = [  sqrt2*[real(H(svones));imag(H(svones)]; diag(H) ]
%         v = [  sqrt2*[H(svones)); diag(H) ]


if strcmp(varflag,'realCase')
    sqrt2=sqrt(2);
    n=length(H);
    %svones: linear indices of **strictly** upper triang. by cols
    svones = logical(triu(ones(n),1));
    v=[sqrt2*H(svones);diag(H)]; % removes 0i
    
elseif strcmp(varflag,'complexCase')
    
    sqrt2=sqrt(2);
    n=length(H);
    %svones: linear indices of **strictly** upper triang. by cols
    svones = logical(triu(ones(n),1));
    v=real([sqrt2*[real(H(svones));imag(H(svones))];diag(H)]); % removes 0i
end
%RH=real(H);
%IH=imag(H);
%v=real([sqrt2*RH(svones);sqrt2*IH(svones);diag(H)]); % removes 0i
end %of function
