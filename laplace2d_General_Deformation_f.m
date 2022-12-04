function [A,B] = laplace2d_General_Deformation_f(F,nx,ny,L,D , eta)
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

%Condition de Pression aux limites
for j=1:ny
    k=(j-1)*nx+1;A3(k,:)=0; Y3(k,:)=0; P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pg;
    k=(j-1)*nx+nx;A3(k,:)=0; Y3(k,:)=0;P3(k,:)=0;P3(k,k) = 1.0; B(2*N + k)=Pd;
    
end
%derivee normale 
for i = 2: nx -1
    k=(ny - 1)*nx +i;
    A3(k,:)=0;Y3(k,:)=0;P3(k,:)=0;
    Y3(k,k) = 1.0; Y3(k , k - nx ) = -1.0;
end
amp = 3;
dec = 2;
%Function inf
for i= 1 : nx
    %positive
    fx = amp*(sin(i)+1) + dec;
    nb = fx + cos(i)*i/5;
    na = amp*cos(i)/(nb);
    d = sqrt(na^2 + nb^2);
    alp = round(na/d);
    blp = round(nb/d);
    
    m = round(fx);
    if m == 0
        m = 1;
    end
    disp(m)
    k = (m-1)*nx + i;
    A1(k,:) = 0;P1(k,:) = 0;A1(k,k) = 1;
    A3(k,:) = 0 ; Y3(k,:)=0;P3(k,:)=0; 
    dirY = 1;
    if blp < 0
        dirY = -1;
    end
    A3(k,k) = alp ;A3(k,k + dirY ) = -alp; % sign plus pour l'autre coté
    Y3(k,k) = blp ; Y3(k,k + nx) = -blp ;
    %derive normale
    
    for j = 1: m-1
        k = (j-1)*nx + i;
        A1(k,:) = 0;P1(k,:) = 0;A1(k,k) = 1;
        A3(k,:) = 0;Y3(k,:) = 0;P3(k,:) = 0;P3(k,k) = 1;
        display (k)
    end

end





Y2 = A1;

A = [A1 Y1 P1; A2 Y2 P2; A3 Y3 P3];