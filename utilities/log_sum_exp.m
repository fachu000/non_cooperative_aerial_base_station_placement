function output = log_sum_exp(input , dim )

% output = log(sum(exp(v_input))) with a trick for numerical stability
% dim : optional parameter. It specifies the dimension along which the sum
%       takes place.

M = max(input,[],'all');
if nargin < 2
	output = log(sum(exp(input-M)))+M;
else
	output = log(sum(exp(input-M),dim))+M;
end

end

