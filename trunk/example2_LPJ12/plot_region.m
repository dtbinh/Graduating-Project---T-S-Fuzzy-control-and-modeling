function plot_region(xi, n, color)

ns = n;
diagonal = zeros(ns,2^ns);
for i = 1:ns
    for k = 0:2^i-1
        for j = k*(2^ns)/(2^i)+1:(k+1)*(2^ns)/(2^i)
            if (k/2 - fix(k/2)) == 0
                diagonal(i,j) = 1;
            else
                diagonal(i,j) = -1;
            end;
        end; % end for j
    end; % end for k
end;  %end for i

for i = 1:2^n
    Chi(:,i) = diag(diagonal(:,i))*xi;
    % diag(v) returns a square diagonal matrix with the elements of vector
    % v on the main diagonal.
end

%plots the polytope that represents the limts of state variables.
line([Chi(1,1);Chi(1,3);Chi(1,4);Chi(1,2);Chi(1,1)],[Chi(2,1);Chi(2,3);...
      Chi(2,4);Chi(2,2);Chi(2,1)], 'Color', color);


end

