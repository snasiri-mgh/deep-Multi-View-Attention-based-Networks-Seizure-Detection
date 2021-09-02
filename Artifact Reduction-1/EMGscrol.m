function [Sigscrol]=EMGscrol(sig,sources,w_x)
D=pinv(w_x)';
N=size(sig);
Sigscrol{1}=sig;
for i=1:N(1)
    D(:,end-i+1:end)=0;
    Sigscrol{i+1}=D*sources+(ones(N(2)-1,1)*mean(sig'))';
end