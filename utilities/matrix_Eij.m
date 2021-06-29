function Eij = matrix_Eij(dim1,dim2,row,col)
%
% Eij = matrix_Eij(dim1,dim2,row,col)
%    DIM1 x DIM2 matrix with all zeros but the element (ROW,COL)
% Eij = matrix_Eij(dim,row,col)
%    DIM x DIM matrix with all zeros but the element (ROW,COL)
%

if nargin ==3
	Eij = zeros(dim1);
	Eij(dim2,row) = 1;
elseif nargin == 4
	Eij = zeros(dim1,dim2);
	Eij(row,col) = 1;
else
	error('invalid syntax'),
end

end
