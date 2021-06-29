
function pos = get_leg_pos( fig )
% pos is a vector with the position [ x y wide height] of the legend of 
% the figure with handle fig. If no arguments are provided, the current
% figure is used

if nargin <1
	fig = gcf;
end

cv = get(fig,'Children');
pos = get(cv(1),'Position');
% 				legpos(3)=1.2*legpos(3);
% 				set(cv(1),'Position',legpos);



end