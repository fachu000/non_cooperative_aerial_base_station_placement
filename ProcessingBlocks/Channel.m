classdef Channel < ProcessingBlock
	%Channel Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		ch_model = 'freeSpace'; % Possible values
		% - 'freeSpace': 
		
		gainConstant = 1; % Gain at unit distance
		distanceConstant = 1; % [distance units]. The square of this 
		%   constant is added to the square distance between the BS and the
		%   MU in the freeSpace model. This prevents the gain and the
		%   gradient from exploding. 
		
	end
	
	methods
		
		% Constructor
		function obj = Channel(varargin)	
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end	
		end
		
	
		function v_gain = gain(obj,m_positionBSs,v_positionMU)
			% m_positionBSs : 3 x num_BS matrix whose n-th column is the
			%     position of the n-th BS
			% v_positionMU : 3 x 1 vector with the position of the MU
			% v_gain : num_BS x 1 vector whose n-th entry is the gain from
			%     the n-th BS to the MU
			
			switch obj.ch_model
				case 'freeSpace'
					v_gain = obj.correctedGainConstant./(sum((m_positionBSs-v_positionMU).^2,1)+obj.distanceConstant^2)'; 
				otherwise
					error('not implemented');
			end
			
		end
		
		function gd = gainGradient(obj,v_positionBS,v_positionMU) 
			% v_positionBS : 3 x 1 vector with the position of the AirBS
			% v_positionMU : 3 x 1 vector with the position of the MU
			% gd : 3 x 1 vector with the gradient of the channel gain.
						
			switch obj.ch_model
				case 'freeSpace'
					gd = 2/obj.correctedGainConstant*(v_positionMU-v_positionBS).*obj.gain(v_positionBS,v_positionMU).^2; % gradient w.r.t. v_positionBS		
				otherwise
					error('not implemented');
			end
			
		end
		
		function cc = correctedGainConstant(obj)
			% When distanceConstant is nonzero, gainConstant needs to be
			% corrected so that it represents the actual gain at unit
			% distance. This corrected constant is the output of this
			% function. 
			
			cc = obj.gainConstant*( 1 + obj.distanceConstant^2 ); 
			
		end
		
		
	end
end

