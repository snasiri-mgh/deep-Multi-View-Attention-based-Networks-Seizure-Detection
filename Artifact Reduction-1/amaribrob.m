function [Xb,Q]=amaribrob(X,n);

[m,N]=size(X);
X=X-(mean(X')'*ones(1,N)); %Removing  mean value
Rxx=(X(:,1:N-1)*X(:,2:N)')/(N-1); %Estimation of sample covariance matrix 
%for the time delay p=1 in order to reduce influence of a white noise.
%
 [Ux,Dx,Vx]=svd(Rxx);
 Dx=diag(Dx);
% n=11;
 if n<m, % under assumption of additive white noise and
        %when the number of sources are known or can a priori estimated
  Dx=Dx-real((mean(Dx(n+1:m))));
  Q= diag(real(sqrt(1./Dx(1:n))))*Ux(:,1:n)';
else    % under assumption of no additive noise and when the 
        % number of sources is unknown
   n=max(find(Dx>1e-99)); %Detection the number of sources
   Q= diag(real(sqrt(1./Dx(1:n))))*Ux(:,1:n)';
end;
Xb=Q*X; % prewhitened data