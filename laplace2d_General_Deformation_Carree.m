function [A,B]=laplace2d_General_Deformation_Carree(F,nx,ny,L,D , eta)

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



P2 = P2 * (beta/2)*(-1/eta);
P1 = P1 * (alpha/2) *(-1/eta);
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

%Condition de Pression
for j=1:ny
    k=(j-1)*nx+1;A3(k,:)=0; Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pg;
    k=(j-1)*nx+nx;A3(k,:)=0; Y3(k,:)=0;P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pd;
end




for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A1(k,k+num)=coeff;
        P3(k,k+num)=coeff;
    end
end

%mettre des valeurs nulles a l'interieur d'obstacle et frontieres
for i=abs:abs+long
    for j=1:ord
        k=(j-1)*nx+i;
        A1(k,:) =0;
        P1(k,:) =0;
        P2(k,:) =0;
        A1(k , k) = 1;
    end
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
    for j=1:ord-1
        k=(j-1)*nx + i;
        A3(k,:) = 0;
        Y3(k,:) = 0;
        P3(k,:) = 0;
        P3(k , k) = 1;
    end
    for j=ny-ord+1:ny
        k=(j-1)*nx + i;
        A3(k,:) = 0;
        Y3(k,:) = 0;
        P3(k,:) = 0;
        P3(k , k) = 1;
    end
end

% conditions a la derive normale
for i=2:abs-1
    k=i; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k , k + nx) = -1.0;
       
    k=(ny - 1)*nx +i;
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k , k - nx ) = -1.0;
end
for i=abs+long+1:nx-1
    k=i; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k ,k + nx) = -1.0;

    k=(ny - 1)*nx +i;
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k ,k - nx ) = -1.0;
end


%obstacle derivee normale
for i=abs:abs+long
    k=(ord - 1)*nx +i;
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k ,k + nx) = -1.0;

    k=(ny - ord - 1)*nx +i;
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k ,k - nx ) = -1.0;
    
end


% k=(ny - 1)*nx +abs;
% Y3(k,k) = -1.0; Y3(k , (ny - 2)*nx + abs ) = 1.0;
% k=(ny - 1)*nx +abs+long;
% Y3(k,k) = -1.0; Y3(k , (ny - 2)*nx + abs+long ) = 1.0;

%derivee normale gauche droite obstocale
for j=2:ord-1
    k=(j-1)*nx+abs; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;

    k=(j-1)*nx+abs+long; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    A3(k,k) = 1.0; A3(k , k + 1) = -1.0;
end

for j = ny - ord + 1: ny - 1
    k=(j-1)*nx+abs; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;

    k=(j-1)*nx+abs+long; 
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    A3(k,k) = 1.0; A3(k , k + 1) = -1.0;    
end

% gauche bas
k = abs;
A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;

P3(k,k-2) =a;
P3(k,k-1)=a;
P3(k,k)=-2*(a+b);
P3(k,k+2*nx)= b;
P3(k,k+nx)=b;

% droite bas
k = abs+long;
A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;
P3(k,k-2) =a;
P3(k,k-1)=a;
P3(k,k)=-2*(a+b);
P3(k,k+2*nx)= b;
P3(k,k+nx)=b;

% gauche haut
k = (ny -1)*nx + abs;
A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;

P3(k,k-2) =a;
P3(k,k-1)=a;
P3(k,k)=-2*(a+b);
P3(k,k-2*nx)= b;
P3(k,k-nx)=b;

% droite haut
k = (ny -1)*nx + abs+long;
A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
A3(k,k) = 1.0; A3(k ,k - 1) = -1.0;
P3(k,k-2) =a;
P3(k,k-1)=a;
P3(k,k)=-2*(a+b);
P3(k,k-2*nx)= b;
P3(k,k-nx)=b;

 
% % Conditions aux conrners
% % obstacle bas
% for k = [abs, (ord -1)*nx +abs]
%         A3(k,:) = 0;
%         Y3(k,:) = 0;
%         P3(k,:) = 0;
% P3(k,k-2) =a;
% P3(k,k-1)=a;
% P3(k,k)=-2*(a+b);
% P3(k,k+2*nx)= b;
% P3(k,k+nx)=b;
% display(k)
% end
% for k = [abs+long ,(ord -1)*nx +abs + long]
%         A3(k,:) = 0;
%         Y3(k,:) = 0;
%         P3(k,:) = 0;
% P3(k,k+2) =a;
% P3(k,k+1)=a;
% P3(k,k)=-2*(a+b);
% P3(k,k+2*nx)= b;
% P3(k,k+nx)=b;
% display(k)
% end
% 
% % obstacle haut
% for k = [(ny-1)*nx + abs  , (ny-ord-1)*nx + abs ]
%         A3(k,:) = 0;
%         Y3(k,:) = 0;
%         P3(k,:) = 0;
% P3(k,k-2) =a;
% P3(k,k-1)=a;
% P3(k,k)=-2*(a+b);
% P3(k,k-2*nx)= b;
% P3(k,k-nx)=b;
% display(k)
% end
% for k=[(ny-1)*nx + abs+long, (ny-ord-1)*nx + abs + long]
%         A3(k,:) = 0;
%         Y3(k,:) = 0;
%         P3(k,:) = 0;
% P3(k,k+2) =a;
% P3(k,k+1)=a;
% P3(k,k)=-2*(a+b);
% P3(k,k-2*nx)= b;
% P3(k,k-nx)=b;
% display(k)
% end


% % Conditions aux conrners
% %obstacle bas
% for k = [abs, (ord -1)*nx +abs]
% A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;A3(k, k)=-1.0; A3(k,k+nx-1)=1.0; Y3(k, k)=1.0; Y3(k,k+nx+1)=-1.0;
% end
% for k = [abs+long ,(ord -1)*nx +abs + long]
% A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;A3(k, k)=1.0; A3(k,k+nx+1)=-1.0; Y3(k, k)=1.0; Y3(k,k+nx+1)=-1.0;
% end
% % obstacle haut
% for k = [(ny-1)*nx + abs  , (ny-ord-1)*nx + abs ]
% A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;A3(k, k)=1.0; A3(k,k-nx-1)=-1.0; Y3(k, k)=-1.0; Y3(k,k-nx-1)=1.0;
% end
% for k=[(ny-1)*nx + abs+long, (ny-ord-1)*nx + abs + long]
% A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;A3(k, k)=1.0; A3(k,k-nx+1)=-1.0; Y3(k, k)=1.0; Y3(k,k-nx+1)=-1.0;
%end



Y2 = A1;

A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];

