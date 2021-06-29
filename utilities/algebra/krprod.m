
function C=krprod( A , B )

% computes the khatri-Rao product, which is defined for matrices A and B with the 
% same number of columns as another matrix C whose n-th column is kron(A(:,n),B(:,n))

n = size(A,2);
assert(size(B,2)==n);
C = zeros( size(A,1)*size(B,1) , n );
for k=1:n
	C(:,k) = kron(A(:,k),B(:,k));
end


end