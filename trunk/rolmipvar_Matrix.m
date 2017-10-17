function poly_M = rolmipvar_Matrix(Mi, label, first_simplex, second_simplex,...
                                    third_simplex, simplexes_degrees)

% This function returns the rolmipvar construction of a matrix dependent of
% alpha, theta, and gamma. Originally, the matrix depends only on one of
%the simplexes, so the other ones hve both degree zero. It is also known
%that the matrix has n vertices, each one of that associated to a different
% component of non-null degree simplex.

degrees_first_simplex = zeros(1, first_simplex);
degrees_second_simplex = zeros(1, second_simplex);
degrees_third_simplex = zeros(1, third_simplex);

if simplexes_degrees(1) == 1
    
    for i = 1:first_simplex
        degrees_first_simplex = zeros(1, first_simplex);
        degrees_first_simplex(i) = 1;
        M{i} = {degrees_first_simplex, degrees_second_simplex,...
                degrees_third_simplex, Mi(:, :, i)};
    end
elseif simplexes_degrees(2) == 1
    for i = 1:second_simplex
        degrees_second_simplex = zeros(1, second_simplex);
        degrees_second_simplex(i) = 1;
        M{i} = {degrees_first_simplex, degrees_second_simplex,...
                degrees_third_simplex, Mi(:, :, i)};
    end
else
    for i = 1:third_simplex
        degrees_third_simplex = zeros(1, third_simplex);
        degrees_third_simplex(i) = 1;
        M{i} = {degrees_first_simplex, degrees_second_simplex,...
                degrees_third_simplex, Mi(:, :, i)};
    end
end
    

poly_M = rolmipvar(M, label, [first_simplex second_simplex third_simplex], simplexes_degrees);
end