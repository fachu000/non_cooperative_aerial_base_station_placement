classdef KmeansAerialBaseStation < AerialBaseStation
	%AerialBaseStation Summary of this class goes here
	%   Detailed explanation goes here
	
	
	
	
% 	properties
% 		v_position = [0 0 0]';
% 		
% 		power ; % Tx. power in natural units.	
% 		stepSize = 1;
% 		
% 		b_fixedHeight = 1; % set to 0 to also adjust height when navigating.
% 		
% 	end
	
	methods
		
		% Constructor
% 		function obj = AerialBaseStation(varargin)	
% 			if nargin > 0
% 				obj =  assignParametersByName(obj,varargin{:});
% 			end	
% 		end
% 		
% 		function v_obj = CloneAtRandomPositions(templateAerialBaseStation,v_power,regionSideLength,b_sameHeight)
% 			%  v_obj is an length(v_power) x 1 vector containing AerialBaseStation objects
% 			%  at different locations, uniformly distributed on a square of
% 			%  side $regionSideLength$. v_power(n) gives the power of the
% 			%  n-th AerialBaseStation
% 			%
% 			%  b_sameHeight : if 1, then the cloned AirBSs have the same
% 			%  height as $templateAerialBaseStation$. 
% 			
% 			v_obj = templateAerialBaseStation.replicate('power',num2cell(v_power),'',[]);
% 			for ind_BS = length(v_power):-1:1
% 				if b_sameHeight
% 					v_obj(ind_BS,1).v_position = [regionSideLength*rand(2,1);templateAerialBaseStation.v_position(3)];
% 				else
% 					error('not implemented');
% 				end
% 			end			
% 			
% 		end
% 		
% 		function m_position = positionAerialBaseStations(v_aerialBaseStations)
% 			% v_aerialBaseStations is an num_BS-length vector; each entry is a
% 			%     AerialBaseStation. 
% 			% m_position  is a 3 x num_BS matrix whose n-th column is the
% 			%     position of the n-th AerialBaseStation
% 			
% 			for ind_BS = length(v_aerialBaseStations):-1:1
% 				m_position(:,ind_BS) = v_aerialBaseStations(ind_BS).v_position;
% 			end
% 			
% 		end
% 		
% 		function v_power = powerAerialBaseStations(v_aerialBaseStations)
% 			% v_aerialBaseStations is an num_BS-length vector; each entry is an
% 			%     AerialBaseStation. 
% 			% v_power  is a num_BS x 1 vector whose n-th entry is the
% 			%     power of the n-th AerialBaseStation
% 			
% 			for ind_BS = length(v_aerialBaseStations):-1:1
% 				v_power(ind_BS,1) = v_aerialBaseStations(ind_BS).power;
% 			end
% 			
% 		end
% 		
		function obj = updatePosition(obj,channel,v_utilityGradient,m_positionsMUs)
			
			% channel : object of class Channel
			% v_utilityGradient : num_MU x 1 vector with the gradient of
			%     the utility w.r.t. to the power received from the present
			%     BS.
			% m_positionsMUs : 3 x num_MU matrix where the n-th column
			%     indicates the position of the m-th MU.
			%			
			
% 			num_MU = size(m_positionsMUs,2); % number of MUs from which we receive a message (miniBatchSize in other files)
% 						
% 			v_stochasticPositionGradient = zeros(3,1);
% 			for  ind_MU = 1:num_MU
% 				% Gradient of the rx power w.r.t. position
% 				v_powerGradient = obj.power * channel.gainGradient(obj.v_position,m_positionsMUs(:,ind_MU));
% 								
% 				v_stochasticPositionGradient = v_stochasticPositionGradient + v_powerGradient*v_utilityGradient(ind_MU);
% 			end
% 			v_stochasticPositionGradient = v_stochasticPositionGradient/num_MU;
% 			
% 			if obj.b_fixedHeight
% 				v_stochasticPositionGradient(3) = 0;
% 			end
           m_positionsClosestMUs = m_positionsMUs(:,v_utilityGradient);

			
	%		obj.v_position = obj.v_position + obj.stepSize*v_stochasticPositionGradient;
	if ~isempty(m_positionsClosestMUs)
	       obj.v_position = mean(m_positionsClosestMUs,2);
	end
		end
		
	end
end

