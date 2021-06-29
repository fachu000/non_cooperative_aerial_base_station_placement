function T = block_toeplitz( blocks_col , blocks_row )

% blocks_col and blocks_row are cell arrays of matrices, which are the
% blocks. similar to the function TOEPLITZ

M = size(blocks_col{1},1 );
N = size(blocks_col{1},2 );
L = length(blocks_col);

T = zeros(M*L,N*L); 

for k=1:L
	E = linshift( eye(L) ,[0 k-1] );
	T = T + kron(E,blocks_row{k});
end

for k=2:L
	E = linshift( eye(L) ,[0 -k+1] );
	T = T + kron(E,blocks_col{k});
end



end