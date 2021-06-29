function X_dft = dft_sample( x , omega )
%
%  X_dft is a length(omega) vector with samples of the DFT of x, understood as a
%  signal x[0],x[1],...,x[length(x)-1].
%
%  X_dft[k] = X( exp( j* (omega(k)) )), with k=1,...,length(omega)
%

x_col = be_column(x);
omega = omega(:);

X_dft = exp( - 1j* omega* (0:length(x_col)-1) ) * x_col;


X_dft = X_dft(:);
if ~is_col_vec(x)
	X_dft = X_dft.';
end
	



end