function str_out = my_sprintf( pattern , varargin )

nargs = length(varargin);

pos = strfind(pattern,'%');
dif = pos(2:end)-pos(1:end-1);
if sum(dif==1)>0
	error('not implemented yet for the %% sequence');
end
if length(pos)~=nargs
	error('number of parameters is incorrect')
end


%%%%% EDIT HERE TO ADD A ESCAPE SEQUENCE %%%%
remove_entries = [];
for k=1:nargs
	switch pattern(pos(k)+1)
		case 'v'
			str_now =  vector_to_str( varargin{k} );	
		otherwise
			continue
	end
	pattern = replace_occurrence( pattern , pos(k) , str_now );
	pos(k+1:end) = pos(k+1:end) + length(str_now)-2;
	remove_entries = [remove_entries k];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

not_to_remove = setxor(1:nargs,remove_entries);
varargin = varargin(not_to_remove);

str_out = sprintf(pattern,varargin{:});

end

function out_cell = remove_kth_entry( in_cell )

end

function new_pattern = replace_occurrence( pattern , pos , str )

new_pattern = [ pattern(1:pos-1) str pattern(pos+2:end) ];

end