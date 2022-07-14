function n=myhist(A,L)
ma=min(A(:));
MA=max(A(:));
A=round((A-ma)*(L-1)/(MA-ma+eps));
x=0:L-1;
n=histc(A,x,1);
n=n';
end