#Implementación Resultado preliminar 1
clear
disp("Matriz original")
A=[0.95 0.04 0 0;0.05 0.85 0 0;0 0.1 1 0;0 0.01 0 1]
N=50000;                        # cantidad de datos generados
[Y X]=ya(N);                 #generacion de datos
disp("Matriz aproximada")
Ap=Y/X                       #Aplicacion de lema 4.2.1

disp("Error interpretativo")
ei=norm(Ap-A,'fro')        #Error interpretativo
disp("Error predictivo")
ep=norm(A*X-Ap*X)         #Error predictivo