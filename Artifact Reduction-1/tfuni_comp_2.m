function U=tfuni2(z);
z1=z(1,:);
z1c=conj(z1);
z2=z(2,:);
z2c=conj(z2);

%Calcul des moments d ordre 2
%-----------------------------------------
m11=mean(z1.*z1c);
m11n=mean(z1.*z1);
m22=mean(z2.*z2c);
m22n=mean(z2.*z2);
m12=mean(z1.*z2c);
m12n=mean(z1.*z2);


%Calcul des moments d ordre 4 (le contraste ne fait intervenir que les cumulants
%ayant autant de puissances complexes que de puissances complexes conjuguees)
%------------------------------------------
m1111=mean(z1.*z1.*z1c.*z1c);
m1112=mean(z1.*z1.*z1c.*z2c);
m1212=mean(z1.*z2.*z1c.*z2c);
m1122=mean(z1.*z1.*z2c.*z2c);
m1222=mean(z1.*z2.*z2c.*z2c);
m2222=mean(z2.*z2.*z2c.*z2c);

%Calcul des cumulants a l ordre 4
%------------------------------------------
c1111 = m1111-abs(m11n)^2-2*m11^2;
c1112 = m1112-m11n*conj(m12n)-2*m11*m12;
c1212 = m1212-2*abs(m12)^2-m11*m22;
c1122 = m1122-m11n*conj(m22n)-2*m12^2;
c1222 = m1222-m12n*conj(m22n)-2*m12*m22;
c2222 = m2222-abs(m22n)^2-2*m22^2;

B=zeros(3,3);
%Calcul des composantes de la matrice symétrique B
%-------------------------------------------------
B(1,1)=c1111+c2222;
B(2,1)= real(c1112-c1222);
B(3,1)= imag(c1112-c1222);
B(2,2)=(1/2)*(c1111+c2222) + 2*c1212 + real(c1122);
B(2,3)=imag(c1122);
B(3,3)=(1/2)*(c1111+c2222)+2*c1212-real(c1122);
B(1,2)=B(2,1);
B(1,3)=B(3,1);
B(2,3)=B(3,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Détermination de la matrice U
[V,D]=eig(B);
[a,p]=max(abs(diag(D)));
w=V(:,p);

alpha=acos(w(1))/2;
if w(3)<w(2)
    phi=atan(w(3)/w(2));
else 
    phi=pi/2 -atan(w(2)/w(3));
end 
c= cos(alpha);
s=sin(alpha);
U=[c s*exp(j*phi); -s*exp(-j*phi) c];
