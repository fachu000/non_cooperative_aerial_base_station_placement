
function [F,d]=linear_system_to_affine_notation( A , b )

% This function computes "F" and "d" such that all the solutions of the
% linear system of equations 
%
%               A * x = b
%
% can be written as
%
%               x =  F*alpha + d
%
% for some alpha and, conversely, any alpha leads to a solution of the
% system.
%
% INPUT
%  A     : is an M x N matrix, where N >= M. It is assumed full row rank
%  b     : is an M x 1 vector
%
% OUTPUT
%  F     : is an N x D matrix, where D is the dimension of the solution
%          subspace, which is given by
%                   D = N - rank(A) 
%          The columns of F are orthonormal
%  d     : is an N x 1 vector in case that there is at least one solution,
%          and it is empty if there is no solution
%

% dummy imput
% M = 30;
% N = 60;
% A = randn(M,N);
% b = randn(M,1);
%A = [A A A];
%%

M = size(A,1);
N = size(A,2);
assert(rank(A)==min(M,N));
r = rank(A);

if r<N
	d = (A'/(A*A'))*b;
elseif (r == N)&&(N==M)
	d = A\b; 
else 
	d = zeros(N,0);
	F = zeros(N,0);
	return
end % any particular solution is OK

[~,F] = orthonormalize(A'); % the columns of F span the complementary subspace 
                        % to the columns of A'

%norm(A*(F*randn(size(F,2),1)+d) - b)



end

