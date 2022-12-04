function [A,B]=laplace2d_General_Deformation_Carree_gen(F,nx,ny,L,D , eta , Pg , Pd , defo )

dx=L/(nx-1); dy=D/(ny-1);
a = 1/dx^2;
b= 1/dy^2;
c = -2/dx^2-2/dy^2;
coeff=[a,b,c,b,a];
alpha = 1/dx;
beta = 1/dy;
num = [  -nx, -1, 0, 1, nx];
N = nx*ny;

A1 = zeros(N,N);
A2 = zeros(N,N);
A3 = zeros(N,N);

Y1 = zeros(N,N);
%Y2 = zeros(N,N);
Y3 = zeros(N,N);

P1 = zeros(N,N);
P2 = zeros(N,N);
P3 = zeros(N,N);

B=reshape(F,3*N,1);


%Pression dans domaine
for i=2:nx-1
    for j=2:ny-1
        if defo(i , j) == 1
            k=(j-1)*nx+i;
            %Derivee de Pression par X
            P1(k , k + [-1 ,0, 1]) = [-1, 0 , 1]; 
            %Derivee de Pression par Y
            P2(k , k + [-nx ,0, nx]) = [-1, 0 , 1];
            %Laplacien de X
            A1(k,k+num)=coeff;
        end
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
    %Ux = 0 aux frontieres haut et en bas
    k = i;
    m=0;
    while 1
        if defo(i,j) == 0
            A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
            Y2(k-nx,:) = 0; P2(k-nx,:)=0; Y2(k-nx,k-nx)= 1.0;
            end
            k=k+nx;
            m=m+1;
        else
            break;
        end
    end
    A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
    
    k=(ny-1)*nx+i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
end




% Conditions aux limites
for j=2:ny-1   
    %Ux = cst en droite et en gauche
    k=(j-1)*nx+1;
    A1(k,:)=0;P1(k,:)=0; 
    A1(k,k)=1.0; A1(k,j*nx)=-1.0; 
    % derivee(Ux) = cst en droite et en gauche
    k =j*nx;
    A1(k,:)=0;P1(k,:)=0; 
    A1(k,(j-1)*nx+2)=1.0; A1(k,(j-1)*nx+1)=-1.0;A1(k,k)=-1.0;A1(k,k-1)=1.0;   
end







end

