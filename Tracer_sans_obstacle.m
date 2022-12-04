clc
clear all
%Variables Global
global M N D L eta Pg Pd
M = 120;N = 40;D = 0.0075;L = 0.0532;eta = 1.79e-5;Pg = 0.082;Pd = 0;
F=zeros(M,3*N);
[A,B]=laplace2d_General_v0(F,M,N,L,D,eta,Pg,Pd);
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
X = (0:L/(M-1):L);Y = (0:D/(N-1):D);
figure(1); surfc( X, Y , U_xy.'); title('Norme de Vitesse en (SI)');xlabel('Longeur X'); ylabel('Hauteur Y'); shading interp; colorbar;
figure(2); surfc( X, Y , Pr.'); title('Pression en (SI)');xlabel('Longeur X'); ylabel('Hauteur Y'); shading flat; colorbar;
%Vitesse theorie
U_t = -(Pg - Pd)*(Y - D).*Y*(1/(2*eta*L));
figure(3);plot(Y,U_t);title('Vitesse theorique pour x donnée');xlabel('Hauteur Y'); ylabel('Vitesse');
x0 = 10; % l'abscisse de calcul
figure(4);plot(Y,U_xy(x0 , :));title('Vitesse calculée pour x donnée');xlabel('Hauteur Y'); ylabel('Vitesse');
%Test de validation
[IsValid , ErrM, ErrM1, ErrM2] = Validate_StokesEq(U_x , U_y , Pr);
disp(['Erreur de validation : ', num2str(IsValid) ]);
%Debit
DQ = zeros(1,M);
for k=1:M
    Q = mean(U_xy(k,:)) * D;
    DQ(k)=Q;
end
figure(4);plot(X,DQ);title('Debit calculée');xlabel('Longeur X'); ylabel('Debit');axis([ 0 L 0 1e-2 ]);
