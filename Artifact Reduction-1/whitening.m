function [W,varb]=whitening(R,n)

m=length(R);

[V,D]=eig(R);
[EigVal_Rxx,index]=sort(diag(D));
VV=V(:,index);
DD=EigVal_Rxx;

varb=0;
if m>n
  for i=1:m-n
      varb=varb+DD(i);
  end
  varb=varb/(m-n);
else
    varb=0;
end


W=[];

for i=m-n+1:m
    W=[W 1/sqrt(DD(i)-varb)*VV(:,i)];
end

W=W';

