function alpha = membership_functions(Mij, N)

%alpha = ones(1,2^N);

i = 0;
for j = 1:2
    for k = 1:2
        for l = 1:2
            for p = 1:2
                i = i + 1;
                %funcoes de pertinencia
                alpha(i)  = Mij(j,1)*Mij(k,2)*Mij(l,3)*Mij(p,4);
            end
        end
    end
end


end

