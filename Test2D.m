clear

M = 100;
N = 100;

C = 5;

dx=1/(M-1); dy=1/(N-1);
X=[0:dx:1]; Y=[0:dy:1];


F=ones(M,N) * C;
[A,B]=laplace2d(F,M,N);

U=A\B;


U1=reshape(U,M,N);


XF = zeros(M,N);

for ind=1:N*N
    [i ,j] = reverse_index(ind , M);
    disp(i);
    disp(j);
    XF(i,j)= U(ind);
end



image(XF*1000)
colorbar


surfc(X,Y,U1); title('Vitesse'); shading interp;

