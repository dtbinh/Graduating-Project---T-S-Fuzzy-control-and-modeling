function x_k = StateVariablesVertices(x1_)
    x2_ = x1_;
    l = 1;
    for i = 1:2
        for j = 1:2
            x_k(1, l) = x1_(i);
            x_k(2, l) = x2_(j);
            l = l + 1;
        end
    end

end

