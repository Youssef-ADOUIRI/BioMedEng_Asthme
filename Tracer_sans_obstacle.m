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
<<<<<<< HEAD
 end
 
 X = (0:L/(M-1):L);
 Y = (0:D/(N-1):D);
 
figure(1); surfc( X, Y , U_xy.'); title('Profil de la vitesse');xlabel('x (en m)'); ylabel('y (en m)'); shading interp; 
h1 = colorbar;
set(get(h1,'label'),'string','Vitesse u (en m/s)');
figure(2); surfc( X, Y , Pr.'); title('Profil de la pression');xlabel('x (en m)'); ylabel('y (en m)'); shading flat;
h2 = colorbar;
set(get(h2,'label'),'string','Pression P (en Pa)');

=======
end
X = (0:L/(M-1):L);Y = (0:D/(N-1):D);
figure(1); surfc( X, Y , U_xy.'); title('Norme de Vitesse en (SI)');xlabel('Longeur X'); ylabel('Hauteur Y'); shading interp; colorbar;
figure(2); surfc( X, Y , Pr.'); title('Pression en (SI)');xlabel('Longeur X'); ylabel('Hauteur Y'); shading flat; colorbar;
>>>>>>> c6d99e365f224855fc6f41c7c4e01a831e67b8e6
%Vitesse theorie
U_t = -(Pg - Pd)*(Y - D).*Y*(1/(2*eta*L));
figure(3);plot(Y,U_t);title('Vitesse theorique en fonction de y');xlabel('y (en m)'); ylabel('Vitesse (en m/s)');
x0 = 10; % l'abscisse de calcul
<<<<<<< HEAD
figure(4);plot(Y,U_xy(x0 , :));title('Vitesse numerique en fonction de y');xlabel('y (en m)'); ylabel('u (en m/s)');

%Test de validation 
=======
figure(4);plot(Y,U_xy(x0 , :));title('Vitesse calculée pour x donnée');xlabel('Hauteur Y'); ylabel('Vitesse');
%Test de validation
>>>>>>> c6d99e365f224855fc6f41c7c4e01a831e67b8e6
[IsValid , ErrM, ErrM1, ErrM2] = Validate_StokesEq(U_x , U_y , Pr);
disp(['Erreur de validation : ', num2str(IsValid) ]);
%Debit
DQ = zeros(1,M);
for k=1:M
    Q = mean(U_xy(k,:)) * D;
    DQ(k)=Q;
end
figure(5);plot(X,DQ);title('Debit calcule');xlabel('Longeur X'); ylabel('Debit');axis([ 0 L 0 1e-2 ]);
