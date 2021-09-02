function W=fas(Rx,estim_A)
[N,P]=size(estim_A);
inv_Rx=inv(Rx);        %calcul de la pseudo inverse de Rx
num=inv_Rx*estim_A;
den=diag(estim_A'*inv_Rx*estim_A);
for p=1:P
	W(p,:)=num(:,p)'./den(p);
end
