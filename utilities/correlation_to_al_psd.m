function psd = correlation_to_al_psd(r,N)
% arbitrary length psd from correlation vector. Use correlation_to_psd if
% possible, since it is more efficient.
%
% r is a L x 1 Hermitian autocorrelation vector. L is an odd number so that
% r( (L-1)/2 ) is real and r( (L-1)/2 - n)  = r( (L-1)/2 + n)'  for all n
% 
% psd has N coefficients
% psd is the normalized psd
% (see signal processing notes)

P = 1/(2*pi)*dft_sample_points( r , N );

n0 = (length(r)-1)/2;
assert( mod(n0,1)==0);
omega = 2*pi*(0:N-1)/N;

if is_col_vec(P)
	psd = exp( 1j*omega*n0 ).'.*P;
else
	psd = exp( 1j*omega*n0 ).*P;
end

assert(norm(imag(psd))<1e-10);

psd = real(psd);

end