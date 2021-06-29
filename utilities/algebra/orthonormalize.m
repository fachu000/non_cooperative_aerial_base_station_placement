function [F,G] = orthonormalize( A )
% The columns of F are an orthonormal basis for the subspace spanned by the
% columns of A.  The columns of G span the orthonormal subspace


[U,~,~] = svd(A);

F = U(:,1:rank(A));
G = U(:,rank(A)+1:end);

end