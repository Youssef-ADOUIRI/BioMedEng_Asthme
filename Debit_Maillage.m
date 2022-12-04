clc
clear all

%Variables Global
global M N D L eta Pg Pd
M = 120;
N = 40;
D = 0.0075;
L = 0.0532;
eta = 1.79e-5;
Pg = 0.082;
Pd = 0;

%Variables Local (obstacle)
abs1 = 50;
ord1 = 4;
long1 = 20;
dx=L/(M-1); dy=D/(N-1);


ite_num = 10;
for k=1:ite_num
    c = 1+z*0.4;
    m = n ;n = round(m*c);
    abs=round(c*abs1);
    ord=round(c*ord1);
    long=round(c*long1);
    F = zeros(m,3*n);
    [A,B]=laplace2d_General_v3(F,M,N,L,D,eta,abs,ord,long);
    U =A\B;
    U_x=reshape(U(1:M*N),M,N);
    U_y=reshape(U(M*N+1:2*M*N),M,N);
    Pr = reshape(U(2*M*N+1:3*M*N),M,N);
    
    %Norme de Vitesse
    U_xy = zeros(M,N);
    for i=1:M
        for j=1:N
            U_xy(i,j)= sqrt(U_x(i,j)^2+U_y(i,j)^2);
        end
    end
    
end