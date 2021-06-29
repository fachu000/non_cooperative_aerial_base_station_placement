function val = my_getfield(si,field)

%   returns si.(field), but field can be of the form sub1.sub2 (nested classes/structs)

ind=strfind(field,'.');
if isempty(ind)
	val = si.(field);
elseif length(ind)==1 
	str1=field(1:ind-1);
	str2=field(ind+1:end);
	ff = si.(str1);
	val = ff.(str2);
else
	error('not implemented');
end

end
