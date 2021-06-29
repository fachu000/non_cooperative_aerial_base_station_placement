function str = vector_to_str( vector )

str = '[';
for k=1:length(vector)
	
	if k == length(vector)
		str = [ str sprintf('%g',vector(k) ) ];
	else
		str = [ str sprintf('%g,',vector(k) ) ];
	end
end
str = [str ']'];

end