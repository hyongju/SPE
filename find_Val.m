



k = 0;
for i = 1:size(A,1)
    if A(i,1) < 0.03 && A(i,2) > 0.97
        k = k+1;
        out(k) = i;
    end
end