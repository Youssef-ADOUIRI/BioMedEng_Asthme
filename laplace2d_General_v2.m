function [A,B]=laplace2d_General_v2(F,nx,ny,L,D , eta)

dx=L/(nx-1); dy=D/(ny-1);
a = 1/dx^2;
b= 1/dy^2;
c = -2/dx^2-2/dy^2;
coeff=[a,b,c,b,a];
alpha = 1/dx;
beta = 1/dy;
num = [  -nx, -1, 0, 1, nx];

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


B=reshape(F,3*N,1);
Pg = 0;
Pd = 1;

%Dans domaine
for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        P1(k , k + [-1 ,0, 1]) = [-1, 0 , 1];
    end
end
for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        P2(k , k + [-nx ,0, nx]) = [-1, 0 , 1];
    end
end

% C.L. sur les frontieres
for i=1:nx
    %Ux = 0 en haut et en bas
    k = i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
    k=(ny-1)*nx+i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
end

% conditions aux limites
for j=2:ny-1
    
    %Ux = cst en droite et en gauche
    k=(j-1)*nx+1;
    A1(k,:)=0; A1(k,k)=1.0; A1(k,j*nx)=-1.0; P1(k,:)=0;
    
    % derivee(Ux) = cst en droite et en gauche
    k =j*nx;
    A1(k,:)=0; A1(k,(j-1)*nx+2)=1.0; A1(k,(j-1)*nx+1)=-1.0;A1(k,k)=-1.0;A1(k,k-1)=1.0; P1(k,:)=0;
    
end

P2 = P2 * beta/2;
P1 = P1 * alpha/2;

%A3 = P1;

%Y3 = P2;

%condition de pression
for j=1:ny
    k=(j-1)*nx+1;A3(k,:)=0;Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pg;
    k=(j-1)*nx+nx;A3(k,:)=0;Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pd;
end
% conditions a la derive normale
for i=2:nx-1
    k=i; 
    Y3(k,k) = 1.0; Y3(k ,k + nx) = -1.0;

    k=(ny - 1)*nx +i;
    Y3(k,k) = -1.0; Y3(k , (ny - 2)*nx +i ) = 1.0;
end







P2 = P2 * -1/eta;
P1 = P1 * -1/eta;

for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A1(k,k+num)=coeff;
        P3(k,k+num)=coeff;
    end
end




Y2 = A1;

%

%A = [full(A1), full((Y1)), full(P1);full(( A2)),full( Y2), full(P2); full(A3), full(Y3 ),full((P3)) ];
A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];
%display(A1);
display(P1);
display(P2);
%display(Y2);
display(A3);
display(P3);
display(Y3);
%display(B);