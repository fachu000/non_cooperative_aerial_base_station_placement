
function b=compute_mse(vec,th_v)
% 
% VEC is M x N, where 
%    M is the number of cases (different parameters/estimators)
%    N is the number of realizations for each case
% TH_V is a M x 1 vector with the theoretical value of the parameter to
% estimate.
% B is an M x 1 vector with the Mean Square Error
N=size(vec,2);
b = mean( abs(vec-th_v*ones(1,N)).^2 ,2);

return

