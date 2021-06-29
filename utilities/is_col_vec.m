function d = is_col_vec( v )
%
% returns 1 if V is  a column vector
%

d = (size(v,1)>=1) && (size(v,2)==1);
return