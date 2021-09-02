function [sn,STDsig]=norma(s)
[n,m]=size(s);

if m<n %Cas vecteur colonne
    s=s';
    sn=(diag(std(s').^-1)*s);
    sn=sn';
    STDsig=diag(std(s'));
else %Cas vecteur ligne
    sn=(diag(std(s').^-1)*s);
    STDsig=diag(std(s'));
end

    
