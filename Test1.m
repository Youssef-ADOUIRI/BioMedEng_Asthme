L = 10 ;
N = 100;
C = 5;

deltaX = L / N;




A = zeros(N+1);
A(1,1) = 1;
A(N+1 , N+1) = 1;
for i = 2:N
    for j = 1:N+1
        if abs(i - j) == 1 
            A(i , j) = 1;
        elseif i == j
            A(i,j)=-2;
        end
        
    end
end

B = ones( N+1 , 1);
B = B*(C*deltaX*deltaX);
B(1,1)=0;
B(N+1 , 1) = 1;
F = linsolve(A,B);

disp(F)

XI = ones(1, N+1);
for index = 1:N+1
    XI(1 , index) = (index-1)*deltaX;
end


XT = transpose(F);
plot(XI , XT);
