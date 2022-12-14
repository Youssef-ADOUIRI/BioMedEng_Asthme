clc
clear all

%Variables Global
global M N D L eta Pg Pd
M = 60;
N = 20;
D = 0.0075;
L = 0.0532;
eta = 1.79e-5;
Pg = 0.082;
Pd = 0;

%Variables Local
dx=L/(M-1); dy=D/(N-1);


ite_num = 20;
Q_arr = zeros(ite_num);
X_arr = zeros(ite_num);
for k=1:ite_num
    c = 1+k*0.2;

    m = round(M*c);
    n = round(N*c);
    
    F = zeros(m,3*n);
    [A,B]=laplace2d_General_v0(F,m,n,L,D,eta,Pg,Pd);
    U =A\B;
    U_x=reshape(U(1:m*n),m,n);
    U_y=reshape(U(m*n+1:2*m*n),m,n);
    Pr = reshape(U(2*m*n+1:3*m*n),m,n);
    
    %Norme de Vitesse
    U_xy = zeros(m,n);
    for i=1:m
        for j=1:n
            U_xy(i,j)= sqrt(U_x(i,j)^2+U_y(i,j)^2);
        end
    end
    
    %Debit
    DQ = zeros(1,m);
    for i=1:m
        DQ(i)=mean(U_xy(i,:)) * D;
    end
    
    Qm = mean(DQ);
    disp([ num2str(k) , ') Débit moyen : ' , num2str(Qm) , ' SI'])
    Q_arr(k) = Qm;
    X_arr(k) = m * n ; 
end
%Ploting
figure(1);plot(X_arr(1:k-1),Q_arr(1:k-1) , 'r--o');
title('Variation de debit en faction de maillage');
xlabel('Taille de Maillage (M*N)'); ylabel('Debit');
%%
% 
% $$e^{\pi i} + 1 = 0$$
% 