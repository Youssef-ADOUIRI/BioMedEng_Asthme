function [A,B]=laplace2d_General_v0(F,nx,ny,L,D,eta,Pg,Pd)

dx=L/(nx-1); dy=D/(ny-1);
a = 1/dx^2;
b = 1/dy^2;
c = -2/dx^2-2/dy^2;
coeff=[b,a,c,a,b];
alpha = 1/dx;
beta = 1/dy;
num = [ -nx, -1, 0, 1, nx];

N=nx*ny;
A1 = sparse(N,N);
A2 = sparse(N,N);
A3 = sparse(N,N);

Y1 = sparse(N,N);
Y3 = sparse(N,N);

P1 = sparse(N,N);
P2 = sparse(N,N);
P3 = sparse(N,N);

B=reshape(F,3*N,1);

%Coeffs de P dans Stokes
for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        P1(k , k + [0, 1]) = [1, -1] * alpha;
        P2(k , k + [0, nx]) = [1, -1] * beta;
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

%div u = 0 centre dans le domaine
for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A3(k , k + [ -1, 1 ]) = [-1, 1]* alpha;
        Y3(k , k + [-nx, nx]) = [-1, 1]* beta;
    end
end

%Laplacien de P a gauche
for j=2:ny-1
    i = 2;
    k = (j-1)*nx + i;
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    P3(k,k+num)=coeff;
end

%condition de pression a g et d
for j=1:ny
    k=(j-1)*nx+1;A3(k,:)=0;Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pg;
    k=(j-1)*nx+nx;A3(k,:)=0;Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pd;
end

% conditions a la derive normale
for i=2:nx-1
    k=(ny-1)*nx + i; 
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k , k - nx) = -1.0;
end

for i=2:nx
    k=i;
    Y3(k,:) = 0; A3(k,:) = 0; P3(k,:) = 0;
    Y3(k,k + [ -nx, -1, 0, 1, nx] + nx)=coeff * eta;
    P3(k,k) = -beta;
    P3(k,k+nx) = beta;
end


for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A1(k,:)=0;
        A1(k,k+num)=eta*coeff;
    end
end


Y2 = A1;
A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];
end
