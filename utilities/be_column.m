function y = be_column( x )

if is_col_vec(x)
	y = x;
elseif is_row_vec(x)
	y = x.';
else
	error(' x is not a vector ');
end

end