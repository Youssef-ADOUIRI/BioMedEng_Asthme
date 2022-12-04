clc 
D=0.00754299;
L=0.0532446375;
Pg = 0.082;
Pd = 0;
eta=1.79e-5;


List_Qt = zeros(1,11); List_Qn = zeros(1,11); List_dx = zeros(1,11);
List_p = zeros(1,11); List_ratio = zeros(1,11);
%D=0.005; L=0.075;

for z = 1:10
m = 300;
n = 20;

c = 1+z*0.4;
M = m ;N = round(n*c);
disp(N)
dx=L/(M-1); dy=D/(N-1);
abs=round(c*50);
ord=round(c*10);
long=round(c*50);
F = zeros(M,3*N);
[A,B]=laplace2d_General_v3(F,M,N,L,D,eta,abs,ord,long);
 U =A\B;
 U_x=reshape(U(1:M*N),M,N);
 U_y=reshape(U(M*N+1:2*M*N),M,N);
 Pr = reshape(U(2*M*N+1:3*M*N),M,N);
 U_xy = zeros(M,N);
 
 
 for i=1:M
    for j=1:N
    U_xy(i,j)= sqrt(U_x(i,j)^2+U_y(i,j)^2);
    end
 end
% Q=0;

arr = zeros( 1 ,M);

for i=1:M
    if (i>=abs && i<= abs +long)
        arr(i) = mean(U_xy(i,1:N-ord))*(N - ord)*D /(N);
    end
    if (i < abs || i > abs +long)
        arr(i) = mean(U_xy(i,:))*D;
    end
end
Q_moy = mean(arr(:));
% figure(1)
% plot(1:M, arr)
% disp(100*(max(arr(:))-min(arr(:)))/Q_moy)
% figure(2)
% surf(U_xy); shading interp; colorbar;

D_moy=(D*(abs-1) + ((N-ord-1)*D/N)*(long) + (M-abs-long)*D)/M ;
%D3_moy=((D^3)*(abs-1) + (((N-ord)*D/N)^3)*(long+1) + (M-abs-long)*(D^3))/M ;
%Qt_moy1 = ((Pg-Pd)/(2*eta*L))*(D_moy^3)/6;
%Qt_moy2 = ((Pg-Pd)/(2*eta*L))*(D3_moy)/6;

Qt1 = (Pg - mean(Pr(abs,1:N-ord-1)))*(D^3)/(12*eta*abs*dx);
Qt2 = (mean(Pr(abs+1,1:N-ord-1)) - mean(Pr(abs+long,1:N-ord-1)))*(D^3)/(12*eta(long-1)*dx);
Qt3 = (mean(Pr(abs+long,1:N-ord-1)) - Pd)*(D^3)/(12*eta(M-abs+long+1)*dx);

Qt_moy = (Qt1+Qt2+Qt3)/3;

disp(Q_moy);
disp(Qt_moy);
disp(Qt_moy2);
p = 100*(max(arr(:))-min(arr(:)))/Q_moy;
ratio = Q_moy/Qt_moy;

List_Qt(1 + z) = Qt_moy; 
List_Qn(1 + z) = Q_moy; 
disp(Qt_moy);
List_dx(1 + z) = dx;
List_ratio(1 + z) = ratio;
List_p(1 + z) = p;
end
disp(List_Qn);
disp(List_dx);
figure(1)
plot(List_dx,List_Qt);
figure(2)
plot(List_dx,List_Qn);
figure(3)
plot(List_dx,List_ratio);
figure(4)
plot(List_dx,List_p);
R1 = 12*eta*((abs-1)*dx/(D^3));
R2 = 12*eta*(long*dx/((N - ord-1)*dy)^3);
R3 = 12*eta*((M-abs-long+1)*dx/(D^3));

 

Rg = R1+R2+R3;
Rn =(Pg-Pd)/Q_moy ;
disp(Rg)
disp(Rn)