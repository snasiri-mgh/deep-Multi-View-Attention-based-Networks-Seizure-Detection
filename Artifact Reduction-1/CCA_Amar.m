
function [S1,W1]=CCA_Amar(sig,Tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(sig);
for i=1:N(1)
sig(i,:)=sig(i,:)-mean(sig(i,:));
end
x=sig(:,1:end-Tau);
y=sig(:,Tau+1:end);
Cxx=(x*x')/(length(x));
Cxy=(x*y')/(length(x));
Cyx=(y*x')/(length(x));
Cyy=(y*y')/(length(x));
[U1,D]=eig(inv(Cxx)*Cxy*inv(Cyy)*Cyx);
[D,INDEX]=sort(diag(D),'descend');
W1=U1';
S1=W1*sig;
S1=S1(INDEX,:);