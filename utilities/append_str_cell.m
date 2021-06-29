
function cello=append_str_cell(cell,str)
% cello=append_str_cell(cell,str)
for k=1:length(cell)
	cstr = cell{k};
	cello{k} = [cstr str];
end
return
