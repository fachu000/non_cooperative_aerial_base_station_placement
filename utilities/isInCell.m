function d=isInCell(str,cellar)
%   d = isInCell(str,cellar)
% returns 0 if STR is not in the cell array of strings CELLAR, and the
% index of the first appearance otherwise.
d=0;
for k=1:length(cellar)
	if strcmpi(str,cellar{k})
		d=k;
		return;
	end
end
end
