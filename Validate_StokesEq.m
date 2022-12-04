function [isValide , ErrM , ErrM1, ErrM2 ] = Validate_StokesEq(X ,Y ,P)
isValide1 = 0;
isValide2 = 0;

global M N D L eta
deltaX = L/(M - 1);
deltaY = D/(N - 1);
ErrM = zeros(M , N); %Matrice d'erreur
ErrM1 = zeros(M,N);
ErrM2 = zeros(M,N);

for i = 2:M-1
    for j = 2:N-1
        %Sur X
        Xx = (X(i+1,j) + X(i-1,j) - 2 * X(i,j))/(deltaX^2);
        Xy = (X(i,j+1) + X(i,j-1) - 2 * X(i,j))/(deltaY^2);
        PX = (P(i+1,j) - P(i,j))/deltaX;
        err1 = eta*(Xx + Xy) - PX;
        ErrM1(i,j) = err1;
        isValide1 = isValide1 + err1; %~=0 si validee
        
        %Sur Y
        Yx = (Y(i+1,j) + Y(i-1,j) - 2 * Y(i,j))/(deltaX*deltaX);
        Yy = (Y(i,j+1) + Y(i,j-1) - 2 * Y(i,j))/(deltaY*deltaY);
        PY = (P(i,j+1) - P(i,j))/deltaY;
        err2 = eta*(Yx + Yy) - PY;
        ErrM2(i,j) = err2;
        isValide2 = isValide2 + err2 ; %~=0 si validee
        
        ErrM(i,j) = sqrt(err1^2 + err2^2);
    end
end

isValide = (isValide1 + isValide2)/((M-1)*(N-1));

end

