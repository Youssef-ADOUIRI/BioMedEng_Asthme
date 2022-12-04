function [A,B]=laplace2d_General_v3(F,nx,ny,L,D , eta, abs , ord , long)


dx=L/(nx-1); dy=D/(ny-1);
a = 1/dx^2;
b= 1/dy^2;
c = -2/dx^2-2/dy^2;
coeff=[a,b,c,b,a];
alpha = 1/dx;
beta = 1/dy;
num = [ -nx, -1, 0, 1, nx];



N=nx*ny;
A1 = zeros(N,N);
A2 = zeros(N,N);
A3 = zeros(N,N);

Y1 = zeros(N,N);
Y3 = zeros(N,N);

P1 = zeros(N,N);
P2 = zeros(N,N);
P3 = zeros(N,N);


B=reshape(F,3*N,1);
Pg = 0.082;
Pd = 0;

%Coeffs de P dans Stokes
for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        P1(k , k + [ 0,1]) = [-1,  1];
        P2(k , k + [0, nx]) = [-1,  1];
    end
end

% C.L. sur les frontieres
for i=1:nx
    %Ux = 0 en haut et en bas
    k = i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
    k=(ny-1)*nx+i; A1(k,:) = 0; P1(k,:)=0; A1(k,k)= 1.0;
end

% % conditions aux limites
for j=2:ny-1
    %Ux = cst en droite et en gauche
    k=(j-1)*nx+1;
    A1(k,:)=0; A1(k,k)=1.0; A1(k,j*nx)=-1.0; P1(k,:)=0;
    
    % derivee(Ux) = cst en droite et en gauche
    k =j*nx;
    A1(k,:)=0; A1(k,(j-1)*nx+2)=1.0; A1(k,(j-1)*nx+1)=-1.0;A1(k,k)=-1.0;A1(k,k-1)=1.0; P1(k,:)=0;
    
end

P2 = P2 * beta;
P1 = P1 * alpha;

%div u = 0 centre dans le domaine
for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A3(k , k + [-1, 1]) = [-1, 1];
        Y3(k , k + [-nx, nx]) = [-1, 1];
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
    %Y3(k,k+[ -nx, -1, 0, 1, nx])=coeff;
    %P3(k,k)=beta/eta;P3(k,k+1)=-alpha/eta;
end

% for i=2:nx-1
%     j = 2;
%     k = (j-1)*nx + i;
%     Y3(k,:)=0;
%     A3(k,:)=0;
%     P3(k,k+num)=coeff;
%     P3(k,k)=-beta/eta;P3(k,k+nx)=beta/eta;
% end



A3 = A3 * alpha/2;
Y3 = Y3 * beta/2;

%condition de pression a g et d
for j=1:ny
    k=(j-1)*nx+1;A3(k,:)=0;Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pg;
    k=(j-1)*nx+nx;A3(k,:)=0;Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pd;
end


% conditions a la derive normale
for i=2:abs-1
    k=(ny-1)*nx + i; 
    Y3(k,k) = 1.0; Y3(k , k - nx) = -1.0;
end

for i=abs+long+1:nx-1
    k=(ny-1)*nx +i; 
    Y3(k,k) = 1.0; Y3(k ,k - nx) = -1.0;
end

%obstacle derivee normale
for i=abs+1:abs+long-1
    k=(ny-ord - 1)*nx +i;
    Y3(k,k) = 1.0; Y3(k ,k - nx) = -1.0;
end

for i=2:nx
k=i;
Y3(k,:) = 0; A3(k,:) = 0; P3(k,:) = 0;
% Y3(k,k) = 1.0; Y3(k ,k + nx) = -1.0;
Y3(k,k+[ -nx, -1, 0, 1, nx]+nx)=coeff;
P3(k,k)=beta/eta;P3(k,k+2*nx)=-beta/eta;
% k = i + nx;
% Y3(k,:) = 0; A3(k,:) = 0; P3(k,:) = 0;
% P3(k,k+num)=coeff;
end






%derivee normale gauche droite obstocale
for j=ny-ord+1:ny-1
    k=(j-1)*nx+abs; 
    A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;

    k=(j-1)*nx+abs+long; 
    A3(k,k) = -1.0; A3(k , k + 1) = 1.0;
end



P2 = P2 * -1/eta;
P1 = P1 * -1/eta;

for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A1(k,k+num)=coeff;
    end
end

%obstacle bas


%mettre des valeurs nulles a l'interieur d'obstacle
for i=abs:abs+long
    for j=ny-ord:ny
        k=(j-1)*nx+i;
        A1(k,:) =0;
        P1(k,:) =0;
        P2(k,:) =0;
        A1(k , k) = 1;
    end
end


%pression nulle a l'interieur d'obstacle
for i=abs+1:abs+long-1
    for j=ny-ord+1:ny
        k=(j-1)*nx + i;
        A3(k,:) = 0;
        Y3(k,:) = 0;
        P3(k,:) = 0;
        P3(k , k) = 1;
    end
end

%droite d'obstacele pression
for j=ny-ord:ny - 1
i = abs + long ;
k = (j-1)*nx + i;
    A3(k,:)=0;
    Y3(k,:)=0;
    P3(k,:)=0;
    %P3(k,k+num)=coeff;
    A3(k,k+[ -nx, -1, 0, 1, nx]+1)=coeff;
    P3(k,k)=alpha/eta;P3(k,k+1)=-alpha/eta;
end


% for i=abs:abs + long
% j =ny-ord-1;
%   k = (j-1)*nx + i;
%     Y3(k,:)=0;
%     %Y3(k,k-nx) = -1;
%     P3(k,k+num)=coeff;
% end

%conditions corners immerges
for k=[(ny - ord - 1)*nx + abs,  (ny - ord-1)*nx+abs + long]
    Y3(k,k) = 1.0; Y3(k , k - nx) = -1.0;
end


%conditions corners parois
for k=[(ny - 1)*nx + abs,  (ny -1)*nx+abs + long]
    Y3(k,k+1) =a;
    Y3(k,k-1)=a;
    Y3(k,k)=c;
    Y3(k,k-2*nx)=b;
    Y3(k,k-nx)=b;
    
    P3(k,k)=beta/eta;
    P3(k,k-nx)=-beta/eta;
end

Y2 = A1;
A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];
