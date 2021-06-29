function B = block_diag( blocks )

M = size(blocks{1},1);
N = size(blocks{1},2);
L = length(blocks);

B = zeros(M*L,N*L);
for k=1:L
	B = B + kron( matrix_Eij(L,L,k,k) , blocks{k} );	
end


end