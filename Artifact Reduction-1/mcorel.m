function R=mcorel(x,y,tau)

% Calcul de la matrice de correlation d'un vecteur au retard tau
% Estimation non biaisee

[nx,mx]=size(x);
[ny,my]=size(y);
if nx>mx ; x=x.'; end % Indice temporel en colonne
if ny>my ; y=y.'; end

N=length(x);

sum=0;

for i=1:N-tau

sum=sum+x(:,i)*y(:,i+tau)';

end
R=1/(N-tau+1)*sum;
