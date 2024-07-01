function C=CMat(c)
%function C=CMat(c)
%INPUT:   c 2n^2 vector real and imaginary parts
%OUTPUT:  C n by n matrix
n=round(sqrt(length(c)/2));
C= zeros(n);
C(:) = c(1:n^2) + 1i*c(n^2+1:end);
end  % end of function
