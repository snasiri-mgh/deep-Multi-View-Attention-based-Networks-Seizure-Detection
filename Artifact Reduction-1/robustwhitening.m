function[Z,WR,DWR]=robustwhitening(X,n,lags,maxit);
% Estimation d'un ensemble de matrice de covariance

[m,N]=size(X);
R=[];
for i=1:length(lags)
    Rii=mcorel(X,X,lags(i));
    Rii=(Rii+Rii.')/2;
    R=[R Rii];
end
%calcul de la SVD de R
[U,SIGMA,V] = svd(R);
%Calcul des Ri
for i=1:length(lags)
    Ri(:,:,i)=U.'*R(:,(((i-1)*m)+1):m*i)*U;
end
%alpha QQ
alpha=ones(1,length(lags));
for round = 1 : maxit + 1,
    if round == maxit + 1,
       fprintf('No convergence after %d steps\n', maxit);
    break;
    end 
alphaold=alpha;
%Calcul de Rbar
    for i=1:length(lags)
        RRbar(:,:,i)=alphaold(1,i)*Ri(:,:,i);
    end
    Rbar=sum(RRbar,3);
    [Uround,Dround]=eig(Rbar);
    Test=sign(diag(Dround))'
    if Test==ones(1,length(Test));
        fprintf('Convergence after %d steps\n', round); 
    break;
    end 
    %calul de delta
    [eigsort,index]=sort(diag(Dround));
    u=Uround(:,index(1));
    for i=1:length(lags)
        delta(1,i)=u.'*Ri(:,:,i)*u;
    end
    delta=delta/norm(delta);
    alpha=alphaold+delta;
end
%calcul de la matrice symmetric definie positive Rbarx
for i=1:length(lags)
    RRbarx(:,:,i)=alpha(1,i)*R(:,(((i-1)*m)+1):m*i);
end
%Calucul de svd de Rbarx
Rbarx=sum(RRbarx,3);
[Us,SIGMASS]=eig(Rbarx);
%[Us,SIGMASS]=svd(Rbarx);
[EigVal_Rxx,index]=sort(diag(SIGMASS));
EigVect_Rxx=Us(:,index);
% Estimation of the variance of the noise
if m>n
    noisevar=mean(EigVal_Rxx(1:m-n)); % Estimation possible if n>m
else
%If m=n, %no estimation possible, assume noiseless environment 
    noisevar=0; 
end
fact=sqrt(EigVal_Rxx(m-n+1:m)-noisevar);
for i=1:n      
    WR(:,i)=(1/fact(i))*EigVect_Rxx(:,m-n+i);
end
WR=WR';
DWR=pinv(WR);
Z=WR*X; % Whitened observation signals



