function v=Hvec(H,svones,sqrt2)
%function v=Hvec(H,svones,sqrt2)
%INPUT:   H n by n Hermitian
%OUTPUT:  v = [  sqrt2*[real(H(svones));imag(H(svones)];diag(H)]
if nargin < 3
        sqrt2=sqrt(2);
	if nargin < 2
                n=length(H);
                %svones: linear indices of **strictly** upper triang. by cols
                svones = logical(triu(ones(n),1));
	end
end
%RH=real(H);
%IH=imag(H);
%v=real([sqrt2*RH(svones);sqrt2*IH(svones);diag(H)]); % removes 0i
v=real([sqrt2*[real(H(svones));imag(H(svones))];diag(H)]); % removes 0i
end %of function
