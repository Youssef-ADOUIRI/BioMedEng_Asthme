clc

%Variables Global
global M N D L eta Pg Pd

%Variables Local 
abs1 = 50;
ord1 = 2;
long1 = 5;
dx=L/(M-1); dy=D/(N-1);

ite_num = 50;
R_arr = zeros(ite_num);
Rt_arr = zeros(ite_num);
X_arr = zeros(ite_num);
for k=1:ite_num
    c = 1+k*0.4;
    abs=abs1;
    ord=round(c*ord1);
    if ord >= (N-4)
        break;
    end
    long=long1;
    F = zeros(M,3*N);
    [A,B]=laplace2d_General_v1(F,M,N,L,D,eta,Pg,Pd,abs,ord,long);
    U =A\B;
    U_x=reshape(U(1:M*N),M,N);
    U_y=reshape(U(M*N+1:2*M*N),M,N);
    Pr = reshape(U(2*M*N+1:3*M*N),M,N);
    
    %Norme de Vitesse
    U_xy = zeros(M,N);
    for i=1:M
        for j=1:N
            U_xy(i,j)= sqrt(U_x(i,j)^2+U_y(i,j)^2);
        end
    end
    %Debit
    DQ = zeros(1,M);
    
    for i=1:M
        if i < abs1-1 || i > abs1 + long1+1
            Q = mean(U_xy(i,:)) * D;
        elseif i > abs1 && i < abs1 + long1
            Q = mean(U_xy(i, 1:N - ord1)) * ( N - ord1)*D/(N - 1);
        end
        DQ(i)=Q;
    end
    %Resistance
    Rh = (Pg - Pd)/mean(DQ);
    disp([ num2str(k) , ') Resistance hydraulique : ' , num2str(Rh) , ' SI'])
    Rt = (12*eta*(L-long1*dx)/(D^3) + 12*eta*long1*dx/(D-ord*dy)^3);
    R_arr(k) = Rh;
    Rt_arr(k) = Rt;
    X_arr(k) = ord*dy;
end

figure(1);plot(X_arr(1:k-1),R_arr(1:k-1), X_arr(1:k-1),Rt_arr(1:k-1), 'r--o');
%plot(,'r--o');
legend('Resistance','Resistance Poiseuille')
title('Resistance hydraulique en fonction de hauteur obstacle');
<<<<<<< HEAD
xlabel('epsilon 2 (en m)'); ylabel('Resistance hydraulique (en Pa.s/m²)');
=======
xlabel('E1 (en m)'); ylabel('Resistance hydraulique (en Pa.s/mÂ²)');
>>>>>>> c6d99e365f224855fc6f41c7c4e01a831e67b8e6
