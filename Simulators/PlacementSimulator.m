classdef PlacementSimulator
	%PLACEMENTSIMULATOR Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		num_steps = 20; % number of time slots where nodes send a message
		%  through the reporting channel and the AirBSs update their
		%  position.
		
		miniBatchSize = 1; % number of measurements per update
		%   if NaN, then the miniBatch equals the entire batch. MUs are
		%   always numbered in the same order. 
		
		v_mobileUsers;
		v_aerialBaseStations; % num_BS x 1 vector of AerialBaseStation objects.
		channel; % object of class Channel
		
		pathPlotter % object of class PathPlotter
		
		b_debug = 0;
		
	end
	
	methods % constructor
		
		function obj = PlacementSimulator(varargin)
			% varargin is a sequence of (property_name,property_value)
			% pairs
			% Example:
			%     car = ProcessingBlock('numberOfWheels',4,'year',2009,'color','blue')
			obj =  assignParametersByName(obj,varargin{:});
			
		end
		
	end
	
	methods % Simulation methods
		
		
		function [t_path,v_F] = generateSamplePath(obj)
			%
			% t_path : 2 x num_BS x num_steps tensor with the spatial
			%     coordinates for each AirBS at each step
			% v_F : vector of objects GFigure
					
			num_BS = length(obj.v_aerialBaseStations);
			num_MU = length(obj.v_mobileUsers);
			t_path = NaN(3,num_BS,obj.num_steps);
			if nargout>1 && isempty(obj.pathPlotter)
				error('A pathPlotter is required to obtain the figures');
			end
			
			% Users that send control messages in each time step
			for ind_step = obj.num_steps:-1:1
				if isnan(obj.miniBatchSize)
					m_ind_selectedMUs(:,ind_step) = 1:num_MU;
				else
					m_ind_selectedMUs(:,ind_step) = randperm( num_MU , obj.miniBatchSize )';
				end
			end
			
			% Figure with power and utility function		
			if ~isempty(obj.pathPlotter)
				obj.pathPlotter.initialize(obj.v_mobileUsers(1));
			end
			
				
			for ind_step = 1:obj.num_steps
				
				m_positionsBSs = obj.v_aerialBaseStations.positionAerialBaseStations();			
				t_path(:,:,ind_step) = m_positionsBSs;
				
				if ~isempty(obj.pathPlotter)
					obj.pathPlotter.receiveNewPlacement(obj.v_aerialBaseStations,...
						obj.channel,obj.v_mobileUsers);
				end
				
				% 3.1. A subset of MUs send messages through the reporting
				% channel
				v_ind_selectedMUs = m_ind_selectedMUs(:,ind_step);
				% Gather information to be communicated through the
				% reporting channel
				v_txPower = obj.v_aerialBaseStations.powerAerialBaseStations();
				for ind_selectedMU = length(v_ind_selectedMUs):-1:1
					
					selectedMU = obj.v_mobileUsers(v_ind_selectedMUs(ind_selectedMU));
					% Utility gradients w.r.t. power from MUs
					v_rxPower = v_txPower.*obj.channel.gain(m_positionsBSs,selectedMU.v_position);
					m_utilityGradient(:,ind_selectedMU) = selectedMU.utilityGradient(v_rxPower);
					%    One row per base station, one column per MU
					
					if obj.b_debug
						m_rxPower(:,ind_selectedMU) = v_rxPower; %one column per MU, one row per BS
					end					
				end
				m_positionsSelectedMUs = obj.v_mobileUsers(v_ind_selectedMUs).positionMobileUsers;
				
				if obj.b_debug
					m_ind_selectedMUs(:,ind_step)
					m_utilityGradient
					m_rxPower_dBm_hard = Watt_to_dBm(max(m_rxPower,[],1));
					m_rxPower_dBm_soft = Watt_to_dBm(log_sum_exp(m_rxPower,1));
					for kk = 1:size(m_rxPower,2)
						m_rxPower_dBm_soft2(kk) = Watt_to_dBm(log_sum_exp(m_rxPower(:,kk)));
					end
					hss2 = [m_rxPower_dBm_hard; m_rxPower_dBm_soft; m_rxPower_dBm_soft2]
				end
				
				% 3.2. The navigation module of each AirBS receives these
				% messages and decides on a new location
				for ind_BS = 1:num_BS
					obj.v_aerialBaseStations(ind_BS).updatePosition( obj.channel,...
						m_utilityGradient(ind_BS,:)',m_positionsSelectedMUs);						
				end
				
				
			end
			
			if ~isempty(obj.pathPlotter)
				
				v_F = obj.pathPlotter.finalize(obj.v_aerialBaseStations,...
						obj.channel,obj.v_mobileUsers);
			end
			
			
			
		end
		
		
	end
end


