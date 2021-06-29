function output = soft_max(v_input)

% output = exp(v_input)/(sum(exp(v_input))) with a trick for numerical stability
%
% soft_max is the gradient of log_sum_exp

M = max(v_input);
output = exp(v_input-M)/(sum(exp(v_input-M)));

end

