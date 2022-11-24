function [Ar errorp errori]=E3recu(N)
  
  Aexact=[0.95 0.04 0 0;0.05 0.85 0 0;0 0.1 1 0;0 0.01 0 1];
#generar los datos
  x0=[1 0 0 0].';  
  X=zeros(4,N);
  Y=zeros(4,N);
  x=x0;
  for i=1:N
    X(:,i)=x;
    Y(:,i)=Aexact*x;
    x=Y(:,i);
  endfor


[U S V]=svd(X);


#Se aplica la parte de la estrategia 1
f=@(A)norm((reshape(A,4,4))*U(:,1:3)*S(1:3,1:3)*(V(:,1:3))'-Y,'fro').^2;

#Se define la restriccion de igualdad
h=@(A)[A(1)+A(2)+A(3)+A(4)-1;A(5)+A(6)+A(7)+A(8)-1;A(9)+A(10)+A(11)+A(12)-1;...
A(13)+A(14)+A(15)+A(16)-1;A(3);A(4);A(9);A(10);A(12);A(13);A(14);A(15)];

#Se define las restricciones de desigualdas
g=@(A)[A(1);A(2);A(3);A(4);A(5);A(6);A(7);A(8);A(9);A(10);A(11);A(12);A(13);A(14)...
;A(15);A(16);A(1)-A(5)-A(9)-A(13);A(6)-A(2)-A(10)-A(14);A(11)-A(3)-A(7)-A(15);A(16)-A(4)-A(8)-A(12);...
A(1)-A(2)-A(3)-A(4);A(6)-A(5)-A(7)-A(8);A(11)-A(9)-A(10)-A(12);A(16)-A(13)-A(14)-A(15)];

#Se aplica sqp para resolver el problema
A0=randn(16,1);

[Ar, obj, info, iter, nf, lambda] =sqp(A0,f,h,g,[],[], 300, 10*e^-15);

for i=1:N-3


[U S V]=svd(X(:,1:3+i));  
#Se aplica la parte de la estrategia 1
f=@(A)norm((reshape(A,4,4))*U*S*(V)'-Y(:,1:3+i),'fro').^2;  
[Ar, obj, info, iter, nf, lambda] =sqp(A0,f,h,g,[],[], 300, 10*e^-15);
errorp=norm(reshape(Ar,4,4)*X-Y,'fro')
errori=norm(reshape(Ar,4,4)-Aexact,'fro')
end


E0=[1 0 0 0].';
E1=E0;
for k=1:300
  E0=[E0 Aexact*E0(:,k)];
  E1=[E1 reshape(Ar,4,4)*E1(:,k)];
end

#grafican orbitas
subplot(211);
plot(E0.');
title('Orbitas del modelo original')

subplot(212);
plot(E1.');
title('Orbitas del modelo identificado')


  
  endfunction