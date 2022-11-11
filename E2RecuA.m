#Recuperar la matriz A de caso preliminar1
clear
# Genera los datos
N=500;
[S E]=ya(N);
 Aexact=[0.95 0.04 0 0;0.05 0.85 0 0;0 0.1 1 0;0 0.01 0 1];
 k=1;
 X=E(:,:);
 Y=S(:,:);
#Se define la funcion objetivo
k=1;
f=@(A)norm((reshape(A,4,4))^k*X-Y,'fro').^2;

#Se define la restriccion de igualdad
h=@(A)[A(1)+A(2)+A(3)+A(4)-1;A(5)+A(6)+A(7)+A(8)-1;A(9)+A(10)+A(11)+A(12)-1;...
A(13)+A(14)+A(15)+A(16)-1];

#Se define las restricciones de desigualdas
g=@(A)[A(1);A(2);A(3);A(4);A(5);A(6);A(7);A(8);A(9);A(10);A(11);A(12);A(13);A(14)...
;A(15);A(16)];

#Se aplica sqp para resolver el problema
ite=99;
cont=0;
tol=1*e^(-13);
A0=zeros(16,1);

[Ar, obj, info, iter, nf, lambda] =sqp(A0,f,h,g,[],[], 300, 10*e^-15);
error=norm(Aexact-reshape(Ar,4,4),'fro');


for k=1:ite-1
   k=k+1;
   c=size(X)(2);
   X=E(:,k:c);
   Y=S(:,k:c);
   
   f=@(A)norm((reshape(A,4,4))^k*X-Y,'fro').^2;
   A0=Ar;
   [Ar, obj, info, iter, nf, lambda] =sqp(A0,f,h,g,[],[], 300, 10*e^-15);
   cont=1+cont;
      error=norm(Aexact-reshape(Ar,4,4),'fro');
      if error<tol
        break
      endif
      
endfor
ite=cont;

Ar=reshape(Ar,4,4);
errori=error
errorp=norm(Ar*E-S,'fro')

#genera 300 datos
E0=[1 0 0 0].';
E1=E0;
for k=1:300
  E0=[E0 Aexact*E0(:,k)];
  E1=[E1 Ar*E1(:,k)];
end

#grafican orbitas
subplot(211);
plot(E0.');
title('Orbitas del modelo original')

subplot(212);
plot(E1.');
title('Orbitas del modelo identificado')