function [sources,A_x,autocor]= EMGsorsep(sig)
disp('CCA')
N=size(sig);
for i=1:N(1)
sig(i,:)=sig(i,:)-mean(sig(i,:));
end
x=sig(:,1:end-1);
y=sig(:,2:end);
[Q_x,R_x]=qr(x',0);
[Q_y,R_y]=qr(y',0);
[U,S,V]=svd(Q_x'*Q_y);
sources=(Q_x*U)';
autocor=diag(S);
w_x=(pinv(R_x)*U);
A_x=pinv(w_x)';

