function H=HMat(v,svind,sqrt2)
%function H=HMat(v,svind,sqrt2)
%Input:  v = [v1;v2;d]  vi in R^(t(n)-n)  d in R^n  (REAL NUMBERS)
% Output:  H  Hermitian
n=round(sqrt(length(v)));
if n^2-length(v)>0
	disp('WARNING!!! catastrophic error; HMat index is wrong')
end
if nargin < 3
        sqrt2=sqrt(2);
	if nargin < 2
                %svind: linear indices of **strictly** upper triang. by cols
		%%%% pass this or use matrix representation???
                svind = find(triu(ones(n),1));  %%?????pass in future???
	end
end
tn = n*(n+1)/2; %number of variables in formulation
tnm = tn-n;
H=zeros(n);
H(svind) = (v(1:tnm) + 1i*v(tnm+1:2*tnm))/sqrt2;
H = (H+H');
H(1:n+1:end) = v(2*tnm+1:end);
H = (H+H')/2;
end % of function
