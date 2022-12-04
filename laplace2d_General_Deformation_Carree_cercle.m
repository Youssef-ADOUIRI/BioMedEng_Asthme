function [A,B] = laplace2d_General_Deformation_Carree_cercle(F,nx,ny,L,D , eta)
dx=L/(nx-1); dy=D/(ny-1);
a = 1/dx^2;
b= 1/dy^2;
c = -2/dx^2-2/dy^2;
coeff=[a,b,c,b,a];
alpha = 1/dx;
beta = 1/dy;
num = [  -nx, -1, 0, 1, nx];

abs=4;
ord=4;
long=4;

N=nx*ny;
%A1 = spalloc(N,N,5*N);
A1 = zeros(N,N);
A2 = zeros(N,N);
A3 = zeros(N,N);

Y1 = zeros(N,N);
%Y2 = zeros(N,N);
Y3 = zeros(N,N);

P1 = zeros(N,N);
P2 = zeros(N,N);
P3 = zeros(N,N);


Pg = 1;
Pd = 0; 
%Pression dans domaine
for i=2:nx-1
    for j=2:ny-1
            k=(j-1)*nx+i;
            %Derivee de Pression par X
            P1(k , k + [-1 ,0, 1]) = [-1, 0 , 1]; 
            %Derivee de Pression par Y
            P2(k , k + [-nx ,0, nx]) = [-1, 0 , 1];
            %Laplacien de X
            A1(k,k+num)=coeff;
    end
end
%Laplacien de Y
Y2 = A1;
%Conservation de masse , div(u)=0
P3 = A1;
%Entrer les coefficients dans les matrices des pressions
P2 = P2 * (beta/2)*(-1/eta);
P1 = P1 * (alpha/2) *(-1/eta);


% C.L. sur les frontieres
for i=1:nx
    %Ux = 0 en haut et en bas
    k = i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
    k=(ny-1)*nx+i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
end
%Condition de Pression aux limites
for j=1:ny
    k=(j-1)*nx+1;A3(k,:)=0; Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pg;
    k=(j-1)*nx+nx;A3(k,:)=0; Y3(k,:)=0;P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pd;
end

abs = 4 ;
radius = 6;

%cercle
for i = abs : abs + 2 * radius
    for j = 1 : radius
        distance = sqrt((i - abs + radius)^2 + ( j )^2 );
        if i ~= abs + radius
        if  distance <= radius
            k=(j-1)*nx+i;
            A1(k,:) =0;
            P1(k,:) =0;
            P2(k,:) =0;
            A1(k , k) = 1;
            A3(k,:) = 0;
            Y3(k,:) = 0;
            P3(k,:) = 0;
            P3(k , k) = 1;
        end
        else
           
        end
    end
end


% conditions a la derive normale
for i=2:abs-1
    k=i; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k , k + nx) = -1.0;
end
for i=abs + 2 * radius+1:nx-1
    k=i; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k ,k + nx) = -1.0;
end

for i = 2 : nx-1 
    k=(ny - 1)*nx +i;
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k ,k - nx ) = -1.0;
end






Y2 = A1;

A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];