function [A,B]=laplace2d_General_v4(F,nx,ny,L,D , eta ,Pg,Pd , abs1 , ord1 , long1 , abs2 , ord2 , long2)


dx=L/(nx-1); dy=D/(ny-1);
a = 1/dx^2;
b = 1/dy^2;
c = -2/dx^2-2/dy^2;
coeff=[b,a,c,a,b];
alpha = 1/dx;
beta = 1/dy;
num = [ -nx, -1, 0, 1, nx];

N=nx*ny;
%initialisation des matrices
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
        A3(k , k + [-1, 1]) = [-1, 1]* alpha/2;
        Y3(k , k + [-nx, nx]) = [-1, 1]* beta/2;
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


% conditions a la derive normale hors
for i=2:abs1-1
    k=(ny-1)*nx + i;
    Y3(k,k) = 1.0; Y3(k , k - nx) = -1.0;
end

for i=abs1+long1+1:nx-1
    k=(ny-1)*nx +i;
    Y3(k,k) = 1.0; Y3(k ,k - nx) = -1.0;
end
%derivee normale dans obstacle 1
for i=abs1+1:abs1+long1-1
    k=(ny-ord1 - 1)*nx +i;
    Y3(k,k) = 1.0; Y3(k ,k - nx) = -1.0;
end

% conditions a la derive normale hors 2
for i=2:abs2-1
    k=i;
    Y3(k,:) = 0; A3(k,:) = 0; P3(k,:) = 0;
    % Y3(k,k) = 1.0; Y3(k ,k + nx) = -1.0;
    Y3(k,k+[ -nx, -1, 0, 1, nx]+nx)=coeff*eta;
    P3(k,k)=beta;P3(k,k+2*nx)=-beta;
end
for i=abs2+long2+1:nx-1
    k=i;
    Y3(k,:) = 0; A3(k,:) = 0; P3(k,:) = 0;
    Y3(k,k+[ -nx, -1, 0, 1, nx]+nx)=coeff*eta;
    P3(k,k)=beta;P3(k,k+2*nx)=-beta;
end
%derivee normale dans obstacle 2
for i=abs2+1:abs2+long2-1
    k=(ord2 - 1)*nx +i;
    Y3(k,:) = 0; A3(k,:) = 0; P3(k,:) = 0;
    Y3(k,k+[ -nx, -1, 0, 1, nx]+nx)=coeff*eta;
    P3(k,k)=beta;P3(k,k+2*nx)=-beta;
end



%derivee normale gauche droite obstocale 1
for j=ny-ord1+1:ny-1
    k=(j-1)*nx+abs1;
    A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;
    
    k=(j-1)*nx+abs1+long1;
    A3(k,k) = -1.0; A3(k , k + 1) = 1.0;
end
%derivee normale gauche droite obstocale 2
for j=2:ord2-1
    k=(j-1)*nx+abs2;
    A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;
    
    k=(j-1)*nx+abs2+long2;
    A3(k,k) = -1.0; A3(k , k + 1) = 1.0;
end

for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A1(k,:)=0;
        A1(k,k+num)=coeff*eta;
    end
end

%obstacle bas


%mettre des valeurs nulles a l'interieur d'obstacle 1
for i=abs1:abs1+long1
    for j=ny-ord1:ny
        k=(j-1)*nx+i;
        A1(k,:) =0;
        P1(k,:) =0;
        P2(k,:) =0;
        A1(k , k) = 1;
    end
end
%mettre des valeurs nulles a l'interieur d'obstacle 2
for i=abs2:abs2+long2
    for j=1:ord2
        k=(j-1)*nx+i;
        A1(k,:) =0;
        P1(k,:) =0;
        P2(k,:) =0;
        A1(k , k) = 1;
    end
end

%pression nulle a l'interieur d'obstacle 1
for i=abs1+1:abs1+long1-1
    for j=ny-ord1+1:ny
        k=(j-1)*nx + i;
        A3(k,:) = 0;
        Y3(k,:) = 0;
        P3(k,:) = 0;
        P3(k , k) = 1;
    end
end
%pression nulle a l'interieur d'obstacle 2
for i=abs2+1:abs2+long2-1
    for j=1:ord2-1
        k=(j-1)*nx + i;
        A3(k,:) = 0;
        Y3(k,:) = 0;
        P3(k,:) = 0;
        P3(k , k) = 1;
    end
end

%droite d'obstacele 1 pression
for j=ny-ord1:ny - 1
    i = abs1 + long1 ;
    k = (j-1)*nx + i;
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    A3(k,k+[ -nx, -1, 0, 1, nx]+1)=coeff*eta;
    P3(k,k)=alpha;P3(k,k+1)=-alpha;
end
%droite d'obstacele 2 pression
for j=2:ord2
    i = abs2 + long2;
    k = (j-1)*nx + i;
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    %P3(k,k+num)=coeff;
    A3(k,k+[ -nx, -1, 0, 1, nx]+1)=coeff;
    P3(k,k)=alpha/eta;P3(k,k+1)=-alpha/eta;
end

%conditions corners immerges
for k=[(ny - ord1 - 1)*nx + abs1,  (ny - ord1-1)*nx+abs1 + long1]
    Y3(k,k) = 1.0; Y3(k , k - nx) = -1.0;
end
%conditions corners immerges 2
for k=[(ord2 - 1)*nx + abs2,  (ord2-1)*nx+abs2 + long2]
    Y3(k,k) = 1.0; Y3(k , k + nx) = -1.0;
end

%conditions corners parois
for k=[(ny - 1)*nx + abs1,  (ny -1)*nx+abs1 + long1]
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    
    Y3(k,k+1) =a;
    Y3(k,k-1)=a;
    Y3(k,k)=c;
    Y3(k,k-2*nx)=b;
    Y3(k,k-nx)=b;
    
    P3(k,k)=beta/eta;
    P3(k,k-nx)=-beta/eta;
end
%conditions corners parois 2
for k=[abs2,  abs2 + long2]
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    
    Y3(k,k+1) =a;
    Y3(k,k-1)=a;
    Y3(k,k)=c;
    Y3(k,k+2*nx)=b;
    Y3(k,k+nx)=b;
    
    P3(k,k)=beta/eta;
    P3(k,k+nx)=-beta/eta;
end

Y2 = A1;
A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];

end
