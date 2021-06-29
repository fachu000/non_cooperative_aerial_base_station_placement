function d = ismat( v )
%
% returns 1 if V is strictly a matrix
%

d = (size(v,2)>1) && (size(v,1)>1);
return