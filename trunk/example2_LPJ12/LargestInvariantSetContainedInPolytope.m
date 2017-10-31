function LMIs = LargestInvariantSetContainedInPolytope(LMIs_, x_k, P)

% LMI (9): [1 b_k'; b_k P(alpha)] >= 0.

    LMIs = LMIs_;
    poly = polytope(x_k');
    [H,K] = double(poly);
    [q, ~] = size(H);
    b_k = [];
    for k = 1:q
        b_k = [b_k (H(k, :)/K(k))'];
    end
    
    for k = 1:length(b_k)
        LMIs = [LMIs, [1 b_k(:,k)' ;b_k(:,k) P] >=0];
    end

end

