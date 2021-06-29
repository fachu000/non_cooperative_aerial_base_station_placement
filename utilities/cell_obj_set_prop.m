function so=cell_obj_set_prop(si,prop,values)
%   SO=CELL_OBJ_SET_prop(SI,prop,VALUES)
%  sets the prop prop of elements of the cell of objects SI to the
%  values in vector VALUES
%  
%  input: 
%     SI    : array of structs
%     prop : string with the name of the prop. Up two nesting levels (eg
%             'prop1' or 'prop1.subprop1'
%     VALUES: scalar/vector with the same number of elements as SI
%  output:
%     SO    : equals SI, but SI(k).prop=values(k) if numel(k)>1 and 
%             SI(k).prop=values otherwise
%
so=si;
if (length(so)>1)&&(length(values)==1)
	for k=1:length(si)
		so{k}=my_setprop(so{k},prop,values);
	end
else
	if (length(so)~=length(values))
		error('numel(si) must equal numel(values)');
	end
	for k=1:length(si)
		so{k}=my_setprop(so{k},prop,values(k));
	end	
end

return

function so=my_setprop(si,prop,value)
so=si;
ind=strfind(prop,'.');
if isempty(ind)
	so.(prop)=value;
else 
	str1=prop(1:ind-1);
	str2=prop(ind+1:end);
	ff=si.(str1);
	ff.(str2)=value;
	so.(str1)=ff;
end

return



