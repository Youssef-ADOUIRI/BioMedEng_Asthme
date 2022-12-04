clc
clear all

%Variables Global
global M N D L eta
M = 120;
N = 40;
D = 0.007543;
L = 0.053244;
eta = 1.79e-5;

% %Variables Local
% abs1 = 18;
% ord1 = 1;
% long1 = 1;

F=zeros(M,3*N);
[A,B]=laplace2d_General_v0(F,M,N,L,D,eta );
%disp(size(B));
%disp(size(A));

 U =A\B;
 U_x=reshape(U(1:M*N),M,N); %Vitesse Ux
 U_y=reshape(U(M*N+1:2*M*N),M,N); %Vitesse Uy
 Pr = reshape(U(2*M*N+1:3*M*N),M,N); %Pression
 
 [IsValid , ErrM, ErrM1, ErrM2] = Validate_StokesEq(U_x , U_y , Pr);
%   display(ErrM1);
%   display(ErrM2);
%display(Xlap);
display(IsValid);
%  disp(U_x);
%  disp(U_y);
%  disp(Pr);