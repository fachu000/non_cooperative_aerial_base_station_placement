function X_dft = dft_sample_points( x , Npoints )

%
%  X_dft is a Npoints vector with uniformly-spaced samples of the DFT of x,
%  understood as a signal x[0],x[1],...,x[length(x)-1].
%
%  X_dft[k] = X( exp( j* 2*pi*k / Npoints )), with k=0,...,Npoints-1
%



omega = 2*pi*(0:Npoints-1)/Npoints;

X_dft = dft_sample( x , omega );








end