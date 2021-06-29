function y=fir_filter(Hd,x)
%
% filters a signal with a FIR filter skipping the transient (uses group
% delay) 
%
%   Hd is a filter object (dfilt), like the output of DESIGN, or LP_FILT or
%   a vector with coefficients. 
%

if (isa(Hd,'dfilt.dffir'))
	[~,N]=max(abs(Hd.Numerator)); % length of transient (group delay)
	y=filter(Hd,[x zeros(1,N-1)]);
	y=y(N:end);	
elseif (isa(Hd,'double'))
	[~,N]=max(abs(Hd)); % length of transient (group delay)
	y=filter(Hd,1,[x zeros(1,N-1)]);
	y=y(N:end);
else
	error('wrong type of Hd');	
end

return