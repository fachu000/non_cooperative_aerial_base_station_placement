function obj = assign_parameters_by_name(obj,varargin)
%
% It allows one to create constructors of the form
%
% function obj = constructorName(obj,field,value,field,value...)
% 
% this function is useful for constructors 
%
if nargin>1
	
	if mod(length(varargin),2)~=0
		error('the number of arguments must be odd');
	end
	
	for k=1:2:length(varargin)-1
		obj.(varargin{k})=varargin{k+1};
	end
	
end
end