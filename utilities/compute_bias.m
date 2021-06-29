
function b=compute_bias(vec,th_v)
% 
% VEC is M x N, where 
%    M is the number of cases (different parameters/estimators)
%    N is the number of realizations for each case
% TH_V is a M x 1 vector with the theoretical value of the parameter to
% estimate.
% B is an M x 1 vector with the bias

b = mean(vec,2)-th_v;

return
