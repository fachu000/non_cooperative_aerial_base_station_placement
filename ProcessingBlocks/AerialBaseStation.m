classdef AerialBaseStation < ProcessingBlock
	%AerialBaseStation Summary of this class goes here
	%   Detailed explanation goes here
	
	
	
	
	properties
		v_position = [0 0 0]';
		
		power ; % Tx. power in natural units.	
		
		b_fixedHeight = 1; % set to 0 to also adjust height when navigating.
		
	end
	
	methods
		
		% Constructor
		function obj = AerialBaseStation(varargin)	
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end	
		end
		
		function v_obj = CloneAtRandomPositions(templateAerialBaseStation,v_power,regionSideLength,b_sameHeight,b_2D)
			%  v_obj is an length(v_power) x 1 vector containing AerialBaseStation objects
			%  at different locations, uniformly distributed on a square of
			%  side $regionSideLength$. v_power(n) gives the power of the
			%  n-th AerialBaseStation
			%
			%  b_sameHeight : if 1, then the cloned AirBSs have the same
			%  height as $templateAerialBaseStation$. 
			%  b_2D: if 0, then the y-coordinate of all generated AirBSs
			%  will be 0. 
			if nargin < 5
				b_2D = 1;
			end
			
			v_obj = templateAerialBaseStation.replicate('power',num2cell(v_power),'',[]);
			for ind_BS = length(v_power):-1:1
				if b_sameHeight
					if b_2D
						v_obj(ind_BS,1).v_position = [regionSideLength*rand(2,1);templateAerialBaseStation.v_position(3)];
					else
						v_obj(ind_BS,1).v_position = [regionSideLength*rand(1,1);0;templateAerialBaseStation.v_position(3)];
					end
				else
					error('not implemented');
				end
			end			
			
		end
		
		function m_position = positionAerialBaseStations(v_aerialBaseStations)
			% v_aerialBaseStations is an num_BS-length vector; each entry is a
			%     AerialBaseStation. 
			% m_position  is a 3 x num_BS matrix whose n-th column is the
			%     position of the n-th AerialBaseStation
			
			for ind_BS = length(v_aerialBaseStations):-1:1
				m_position(:,ind_BS) = v_aerialBaseStations(ind_BS).v_position;
			end
			
		end
		
		function v_power = powerAerialBaseStations(v_aerialBaseStations)
			% v_aerialBaseStations is an num_BS-length vector; each entry is an
			%     AerialBaseStation. 
			% v_power  is a num_BS x 1 vector whose n-th entry is the
			%     power of the n-th AerialBaseStation
			
			for ind_BS = length(v_aerialBaseStations):-1:1
				v_power(ind_BS,1) = v_aerialBaseStations(ind_BS).power;
			end
			
		end
		
	end
	
	methods(Abstract)
		
		obj = updatePosition(obj,channel,v_utilityGradient,m_positionsMUs)
		% channel : object of class Channel
		% v_utilityGradient : num_MUFeedback x 1 vector. The n-th entry is the
		%     control message sent by the n-th MU. The total number of
		%     MUs that send a control message may be less than the
		%     number of MUs in the network.
		% m_positionsMUs : 3 x num_MUFeedback matrix where the n-th column
		%     indicates the position of the m-th MU.
		%
		
		
		
	end
end

