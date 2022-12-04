clc 

M = 120;
N = 40;
D=0.007543;
L=0.053244;
eta=1.79e-5;

abs1=50;
ord1=10;
long1=20;

% abs2=30;
% ord2=3;
% long2=15;

Pg = 0.082;
Pd = 0;

F=zeros(M,3*N);
[A,B]=laplace2d_General_v3(F,M,N,L,D,eta, abs1 , ord1 , long1 );
disp(size(B));
disp(size(A));
%display(A);
Z = zeros(3*M*N);
% t=0;
% nb=0;
% for i=1:3*N*M
%     for j=1:3*N*M
%         if A(j,i)~= 0
%             t=1;
%         end
%     end
%     if t==0
%         disp([i,j]);
%         nb = nb+1;
%     end
% t=0;
% end
% disp(nb);

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
Q=0;

arr = zeros( 1, M);
for i=1:N
 U_moy = mean(U_xy(i,:));
 U_moy2 =mean(U_xy(i,1:N-ord1));
 surf = pi*(D/2)^2; % surface hors obstacel
 Q=U_moy*surf;
 if (i>=abs1 && i<= abs1 +long1) 
    surf2 = pi*((D)^2 - (ord1 - D)^2)/4;
    Q = U_moy2*surf2;
 end
disp(Q)
arr(i) = Q;
end
display(arr);
figure(1)
plot(1:M , arr)

figure(3);
surfc((0:D/(N-1):D), (0:L/(M-1):L) ,Pr); title('Pression'); shading interp; colorbar;
figure(4)
surfc( (0:D/(N-1):D), (0:L/(M-1):L) , U_xy); title('Norme de Vitesse'); shading interp; colorbar;

%figure(5)
%surfc(U_y); title('Vitesse uy'); shading interp; colorbar;