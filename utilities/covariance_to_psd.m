function [psd,waxis] = covariance_to_psd(C)
% extract the psd from a toeplitz autocorrelation matrix

c1 = C(:,1).';
c2 = fliplr(C(1,2:end));

r = [ c1 c2 ];
psd = abs(fft( r ));
nl = length(psd);
waxis = (0:nl-1)*(2*pi/nl);



end