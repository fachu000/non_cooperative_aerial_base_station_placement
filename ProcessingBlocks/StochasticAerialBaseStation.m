classdef StochasticAerialBaseStation < AerialBaseStation
	%STOCHASTICAERIALBASESTATION Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		stepSize = 1;
		b_debug = 0;
		
	end
	
	methods
		
		% Constructor
		function obj = StochasticAerialBaseStation(varargin)	
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end	
		end
		
		
		function obj = updatePosition(obj,channel,v_utilityGradient,m_positionsMUs)
			
			% channel : object of class Channel
			% v_utilityGradient : num_MU x 1 vector with the gradient of
			%     the utility w.r.t. the power received from the present
			%     BS.
			% m_positionsMUs : 3 x num_MU matrix where the n-th column
			%     indicates the position of the m-th MU.
			%			
			
			num_MU = size(m_positionsMUs,2); % number of MUs from which we receive a message (miniBatchSize in other files)
						
			v_stochasticPositionGradient = zeros(3,1);
			for  ind_MU = 1:num_MU
				% Gradient of the rx power w.r.t. position
				v_powerGradient = obj.power * channel.gainGradient(obj.v_position,m_positionsMUs(:,ind_MU));
								
				v_stochasticPositionGradient = v_stochasticPositionGradient + v_powerGradient*v_utilityGradient(ind_MU);
			end
			v_stochasticPositionGradient = v_stochasticPositionGradient/num_MU;
			if obj.b_debug
				step = obj.stepSize*v_stochasticPositionGradient
				if norm(step)>100
					keyboard
				end
			end
			if obj.b_fixedHeight
				v_stochasticPositionGradient(3) = 0;
			end
			
% 			ind_MU
% 			v_stochasticPositionGradient
% 			pause
			
			obj.v_position = obj.v_position + obj.stepSize*v_stochasticPositionGradient;
			
		end
		
	end
end

