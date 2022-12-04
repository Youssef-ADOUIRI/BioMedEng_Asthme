clc
clear all

%Variables Global
global M N D L eta Pg Pd
M = 240;
N = 80;
D = 0.0075;
L = 0.0532;
eta = 1.79e-5;
Pg = 0.082;
Pd = 0;

%Variables Local
abs1 = 100;
ord1 = 10;
long1 = 50;
abs2 = 60;
ord2 = 10;
long2 = 20;

F=zeros(M,3*N);
[A,B]=laplace2d_General_v4(F,M,N,L,D,eta,Pg,Pd , abs1 , ord1 , long1, abs2 , ord2 , long2);
disp(size(B));
disp(size(A));

U =A\B;
U_x=reshape(U(1:M*N),M,N); %Vitesse Ux
U_y=reshape(U(M*N+1:2*M*N),M,N); %Vitesse Uy
Pr = reshape(U(2*M*N+1:3*M*N),M,N); %Pression

%Norme de Vitesse
U_xy = zeros(M,N);
for i=1:M
    for j=1:N
        U_xy(i,j)= sqrt(U_x(i,j)^2+U_y(i,j)^2);
    end
end

X = (0:L/(M-1):L);
Y = (0:D/(N-1):D);

figure(1); surfc( X, Y , U_xy.'); title('Profile de vitess');xlabel('x (en m)'); ylabel('y (en m)'); shading interp; h=colorbar;ylabel(h,'Vitesse en m/s')
figure(2); surfc( X, Y , Pr.'); title('Pression');xlabel('Longeur X'); ylabel('Hauteur Y'); shading flat; colorbar;