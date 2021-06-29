
function [A,b]=affine_notation_to_linear_system( F , d )

% This function computes "A" and "b" such that all the solutions of the
% linear system of equations 
%
%               A * x = b
%
% can be written as
%
%               x =  F*alpha + d
%
% for some alpha and, conversely, any alpha below leads to a solution of the
% system above
%
% INPUT
%  F     : is an N x D matrix with independent columns, where D <= N
%  d     : is an N x 1 vector
%
% OUTPUT
%  A     : is an (N-D) x N matrix with orthonormal rows
%  b     : is an (N-D) x 1 vector
%

% dummy imput
% M = 30;
% N = 60;
%%

N = size(F,1);
D = size(F,2);
assert(D<=N);
assert(rank(F)==D);

[~,At] = orthonormalize(F); 
A = At';
b = A*d;

%norm(A*(F*randn(size(F,2),1)+d) - b)



end

