function d = is_complex( A )
%
% returns 1 if A is complex and 0 if it is real
%
thresh = 1e-15;
d = norm(imag(A)) > numel(A)*thresh;


return