clc

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
ord1 = 8;
long1 = 30;

F=zeros(M,3*N);
[A,B]=laplace2d_General_v1(F,M,N,L,D,eta,Pg,Pd , abs1 , ord1 , long1);
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

figure(1); surfc( X, Y , U_xy.'); title('Norme de Vitesse');xlabel('Longeur X'); ylabel('Hauteur Y'); shading interp; colorbar;
figure(2); surfc( X, Y , Pr.'); title('Pression');xlabel('Longeur X'); ylabel('Hauteur Y'); shading flat; colorbar;

%Test de validation
[IsValid , ErrM, ErrM1, ErrM2] = Validate_StokesEq_obs(U_x , U_y , Pr , abs1 , ord1 , long1);
disp(['Erreur de validation : ', num2str(IsValid) ]);

%Debit
DQ = zeros(1,M);
DQt = zeros(1,M);

for k=1:M
    if k < abs1-1 || k > abs1 + long1+1
        Q = mean(U_xy(k,:)) * D;
    elseif k > abs1 && k < abs1 + long1
        Q = mean(U_xy(k, 1:N - ord1)) * ( N - ord1)*D/(N - 1);
    end
    DQ(k)=Q;
end
for k=2:M-1
    if k < abs1-1 || k > abs1 + long1 + 1
        Qt = (mean(Pr(k-1,:)) -  mean(Pr(k+1,:))) * (D^3) / (12*eta*2*L/(M-1));
    elseif k > abs1 && k < abs1 + long1
        d = ( N - ord1 - 1)*D/(N - 1);
        Qt = (mean(Pr(k-1,1:N - ord1 -1)) -  mean(Pr(k+1,1:N - ord1 -1))) * (d^3) / (12*eta*2*L/(M-1));
    end
    DQt(k)=Qt;
end
%fin limite
DQt(1) = (mean(Pr(1,:)) -  mean(Pr(2,:))) * (D^3) / (12*eta*L/(M-1));
DQt(M) = (mean(Pr(M-1,:)) -  mean(Pr(M,:))) * (D^3) / (12*eta*L/(M-1));
p = (max(DQ)-min(DQ))/mean(DQ);
disp(['Incértitude : ' , num2str(p*100) , '%'])
figure(4);plot(X,DQ);title('Debit calculée');xlabel('Longeur X'); ylabel('Debit');axis([ 0 L 0 2.5e-3 ]);
figure(5);plot(X,DQt);title('Debit Thorique');xlabel('Longeur X'); ylabel('Debit');axis([ 0 L 0 3.5e-3 ]);
figure(6);plot(X,(DQt./DQ)*100);title('(Debit Thorique) / (Debit calculée)  ');xlabel('Longeur X'); ylabel('Rapport %');axis([ 0 L 0 100 ]);

%Moyen debit
disp(['Moyen du debit : ' , num2str(mean(Q)) , ' SI'])

%Resistance
Rh = (Pg - Pd)/mean(Q);
disp(['Resistance hydraulique : ' , num2str(Rh) , ' SI'])