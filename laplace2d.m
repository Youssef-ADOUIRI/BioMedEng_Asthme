function [A,B]=laplace2d(F,nx,ny)

dx=1/(nx-1); dy=1/(ny-1);

coeff=[-1/dy^2,-1/dx^2,2/dx^2+2/dy^2,-1/dx^2,-1/dy^2];

num = [  -nx, -1, 0, 1, nx];

N=nx*ny;
A=spalloc(N,N,5*N);
B=reshape(F,N,1);

for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        A(k,k+num)=coeff;
    end
end

% conditions aux limites

for j=1:ny
    k=(j-1)*nx+1; A(k,:)=0; A(:,k);A(k,k)=1.0; B(k)=0;
    k=(j-1)*nx+nx; A(k,:)=0; A(:,k)=0; A(k,k)=1.0; B(k)=0;
end;
% C.L. sur les frontieres
for i=1:nx
    k=i; A(k,:)=0; A(:,k)=0; A(k,k)=1.0; B(k)=0;
    k=(ny-1)*nx+i; A(k,:)=0; A(:,k)=0; A(k,k)=1.0; B(k)=0;
end;

