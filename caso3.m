#Solución del caso 3
function [Ai1 Ai2 B]=caso3(A1,A2,B,N)
 
  X1=randn(5,3*N);
  X2=randn(8,3*N);
  
  Aa=[A1 A2];
  X=[X1;X2];
  Y=Aa*X+repmat(B,1,N);
  Xa=[X;repmat(eye(3),1,N)];
  
  AB=Y/Xa;
  
  Ai1=AB(:,1:5);
  Ai2=AB(:,6:13);
  B=AB(:,14:16);