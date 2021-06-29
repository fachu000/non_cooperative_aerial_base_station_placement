function d = is_row_vec( v )
%
% returns 1 if V is a row vector
%

d = (size(v,2)>=1) && (size(v,1)==1);
return