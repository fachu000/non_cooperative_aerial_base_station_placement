function so=vec_struct_set_field(si,field,values)
%   SO=VEC_STRUCT_SET_FIELD(SI,FIELD,VALUES)
%  sets the field FIELD of elements of the vector of structures SI to the
%  values in vector VALUES
%  
%  input: 
%     SI    : array of structs
%     FIELD : string with the name of the field. Up to two nesting levels (eg
%             'field1' or 'field1.subfield1'
%     VALUES: scalar/vector with the same number of elements as SI
%  output:
%     SO    : equals SI, but SI(k).field=values(k) if numel(k)>1 and 
%             SI(k).field=values otherwise
%

so=si;
if (length(so)>1)&&(length(values)==1)
	for k=1:length(si)
		so(k)=my_setfield(so(k),field,values);
	end
else
	if (length(so)~=length(values))
		error('numel(si) must equal numel(values)');
	end
	for k=1:length(si)
		so(k)=my_setfield(so(k),field,values(k));
	end	
end

return

function so=my_setfield(si,field,value)
ind=strfind(field,'.');
if isempty(ind)
	if ~isfield(si,field)
		error('SI must have a field named FIELD');
	end
	so=setfield(si,field,value);
else 
	str1=field(1:ind-1);
	str2=field(ind+1:end);
	ff=getfield(si,str1);
	ff=setfield(ff,str2,value);
	so=setfield(si,str1,ff);
end

return



