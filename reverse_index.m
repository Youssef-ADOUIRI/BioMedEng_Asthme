function [i,j] = reverse_index(k , N)
i = mod(k,N);
if i == 0
    i = N;
end
j = 1+ (k - i)/N;
end

