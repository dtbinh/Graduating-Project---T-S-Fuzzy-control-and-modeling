function [LMIs, crit] = EnlargementOfLargestInvariantSet(LMIs_, P)

% LMI (9)
LMIs = LMIs_;
beta = sdpvar(1);        
LMIs = [LMIs, beta >=0];
In = eye(2);
LMIs = [LMIs, P - beta*In <=0];
crit = beta;

end

