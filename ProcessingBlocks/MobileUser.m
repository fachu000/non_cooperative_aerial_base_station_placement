classdef MobileUser < ProcessingBlock
	%MOBILEUSER Summary of this class goes here
	%   Detailed explanation goes here
	
% 	properties(Constant)
% 		f_sigmoid = @(x) exp(x)./(exp(x)+1);
% 		f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
% 	end
	
	properties
		v_position = [0 0 0];
		
		ch_utility = ''; % Possible values
		%   - 'sigmoid' --> utility is the shifted sigmoid 
		%   - 'leakyRelu': min( positiveScaling*(p - powerStart) , p - powerStart )
		%      where p is result of the aggregate function
		%   - 'kmeans'
		
		ch_aggregate = 'maxPower';
		%   - 'totalPower' : p is the sum of the power rx from all BSs
		%   - 'smoothMaxPower'   : p is log_sum_exp of the power rx from all BSs
		%   - 'maxPower'   : p is the maximum of the power rx from the BSs
		
		% Utility parameters:
		powerStart  % in Watt, used by the shifted sigmoid and leaky relu
		powerStop   % in Watt, used by the shifted sigmoid
		%   The shifted sigmoid is roughly 0 for inputs less than
		%   powerStart, and roughly 1 for power greater than powerStop.
		
		b_dBm = 0; % If 1, then the class operates internally in dBm, despite
		%   the fact that all input arguments it takes are always in Watt. The
		%   log-sum-exp function for computing maxPower operates in dBm.
		%   The parameters powerStart and powerStop are converted to dBm
		%   and the utilities enforce these limits in dBm.
		
		% Parameter for the leaky relu
		positiveScaling = .2; 
	end
	
	methods
		
		% Constructor
		function obj = MobileUser(varargin)
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end
		end
		
		function v_obj = cloneAtRandomPositions(templateMobileUser,num_MU,regionSideLength,b_sameHeight,b_2D)
			%  v_obj is an num_MU x 1 vector containing MobileUser objects
			%  at different locations, uniformly distributed on a square of
			%  side $regionSideLength$
			%
			%  b_sameHeight : if 1, then the cloned MUs have the same
			%  height as obj.
			%
			%  b_2D: if 0, then all resulting MUs will have 0 y-coordinate
			if nargin < 5
				b_2D = 1;
			end
			
			for ind_MU = num_MU:-1:1
				v_obj(ind_MU,1) = templateMobileUser.clone();
				if b_sameHeight
					if b_2D
						v_obj(ind_MU,1).v_position = [regionSideLength*rand(2,1);templateMobileUser.v_position(3)];
					else
						v_obj(ind_MU,1).v_position = [regionSideLength*rand(1,1);0;templateMobileUser.v_position(3)];
					end
				else
					error('not implemented');
				end
			end
			
		end
		
		function v_obj = cloneAtPositions(templateMobileUser,m_positions)
			%  v_obj is an size(m_positions,2) x 1 vector whose n-th entry
			%  is an object of class MobileUser with v_position
			%  m_positions(:,n)
			%
			 
			num_MU = size(m_positions,2);
			for ind_MU = num_MU:-1:1
				v_obj(ind_MU,1) = templateMobileUser.clone();				
				v_obj(ind_MU,1).v_position = m_positions(:,ind_MU);				
			end
			
		end
		
		function m_position = positionMobileUsers(v_mobileUsers)
			% v_mobileUser is an num_MU-length vector; each entry is a
			%     MobileUser. 
			% m_position  is a 3 x num_MU matrix whose n-th column is the
			%     position of the n-th MobileUser
			
			for ind_MU = length(v_mobileUsers):-1:1
				m_position(:,ind_MU) = v_mobileUsers(ind_MU).v_position;
			end
			
		end
		
		function u = utility(obj,v_rxPower)
			% v_rxPower is a num_BS x 1 vector with the power received from
			%    each BS in Watt.
			
			p = aggregate(obj,v_rxPower);	% p is in dBm if b_dBm == 1, in Watt otherwise
			
			switch obj.ch_utility
				case 'sigmoid'
					u = obj.shiftedSigmoid(p);
				case 'leakyRelu'
					u = obj.shiftedLeakyRelu(p); 
				case 'kmeans'					
					u = NaN; %obj.shiftedSigmoid(p);					
				otherwise
					error('not implemented');
			end
			
% 			if isnan(sum(v_rxPower)) || isnan(u)
% 				keyboard
% 			end
				
			
		end
		
		function p = aggregate(obj,v_rxPower)
			% v_rxPower is a num_BS x 1 vector with the power received from
			%    each BS in Watt.
			switch obj.ch_aggregate
				case 'totalPower'	
					p = sum(v_rxPower);
					if obj.b_dBm
						p = Watt_to_dBm(p);						
					end					
				case 'smoothMaxPower'
					if obj.b_dBm
						v_rxPower = Watt_to_dBm(v_rxPower);						
					end		
					p = log_sum_exp(v_rxPower);
				case 'maxPower'
					if obj.b_dBm
						v_rxPower = Watt_to_dBm(v_rxPower);						
					end	
					p = max(v_rxPower);					
				otherwise
					error('not implemented')
			end
			
		end
			
		
		function v_grad = aggregateGradient(obj,v_rxPower)
			% v_rxPower is a num_BS x 1 vector with the power received from
			%    each BS in Watt.
			% v_grad is a num_BS x 1 vector with the gradient of the
			%    aggregate function w.r.t. the rx power from each BS.
			%
			num_BSs = length(v_rxPower);
			
			switch obj.ch_aggregate
				case 'totalPower'	
					v_grad = ones(num_BSs,1);
					if obj.b_dBm
						p = sum(v_rxPower);
						v_grad = v_grad * Watt_to_dBm_grad(p);						
					end					
				case 'smoothMaxPower'
					v_grad = soft_max(v_rxPower);
					if obj.b_dBm
						v_grad = Watt_to_dBm_grad(v_rxPower).*v_grad; 
					end
				case 'maxPower'
					error('not implemented')
				otherwise
					error('not implemented')
			end
			
		end

		
		function v_gradient = utilityGradient(obj,v_rxPower)
			% v_rxPower is a num_BS x 1 vector with the power received from
			%    each BS in Watt.
			% v_gradient is a num_BS x 1 vector whose m-th entry is the
			%    derivative of the utility w.r.t. the power received from the
			%    m-th BS
			
			p = aggregate(obj,v_rxPower);	% p is in dBm if b_dBm == 1, in Watt otherwise
			v_aggregateGradient = aggregateGradient(obj,v_rxPower);
			
			switch obj.ch_utility
				case 'sigmoid'
					up = obj.shiftedSigmoidPrime(p);					
					v_gradient = v_aggregateGradient*up;
				case 'leakyRelu'
					up = obj.shiftedLeakyReluPrime(p); 
					v_gradient = v_aggregateGradient*up;
				case 'kmeans'					
					up = NaN; %obj.shiftedSigmoid(p);					
					v_gradient = (v_rxPower == max(v_rxPower));
				otherwise
					error('not implemented');
			end
			
		end
		
			
		function gradient = utilityGradientOneBS(obj,v_rxPower,ind_BS)
			% gradient is the derivative of the utility w.r.t. the power received from the
			%    ind_BS-th BS
			
			v_gradient = obj.utilityGradient(v_rxPower);
			gradient = v_gradient(ind_BS);
		end
		
		
		function s = shiftedSigmoid( obj , x )
			% x : power in watt if obj.b_dBm == 0 or in dBm	 otherwise
			if isempty(obj.powerStop)
				error('You must specify a value for obj.powerStop');
			end
			if obj.b_dBm				
				powerStart_h = Watt_to_dBm(obj.powerStart);
				powerStop_h = Watt_to_dBm(obj.powerStop);
			else
				powerStart_h = obj.powerStart;
				powerStop_h = obj.powerStop;
			end
			s = sigmoid(3*(2*x-powerStart_h-powerStop_h)/(powerStop_h-powerStart_h));
			
		end
		
		function s = shiftedSigmoidPrime( obj , x )
			% x : power in watt if obj.b_dBm == 0 or in dBm	 otherwise
			if obj.b_dBm
				powerStart_h = Watt_to_dBm(obj.powerStart);
				powerStop_h = Watt_to_dBm(obj.powerStop);
			else
				powerStart_h = obj.powerStart;
				powerStop_h = obj.powerStop;
			end
			s = sigmoid_prime(3*(2*x-powerStart_h-powerStop_h)/(powerStop_h-powerStart_h))...
				* 3*2/(powerStop_h-powerStart_h);			
		end
		
		function s = shiftedLeakyRelu( obj , x )
			% x : power in watt if obj.b_dBm == 0 or in dBm	 otherwise
			
			if obj.b_dBm					
				powerStart_h = Watt_to_dBm(obj.powerStart);
			else
				powerStart_h = obj.powerStart;
			end
						
			x_shifted = x - powerStart_h;
			s = min( x_shifted ,obj.positiveScaling*x_shifted);
		end
		
		function s = shiftedLeakyReluPrime( obj , x )
			% x : power in watt if obj.b_dBm == 0 or in dBm	 otherwise
			
			if obj.b_dBm
				powerStart_h = Watt_to_dBm(obj.powerStart);
			else
				powerStart_h = obj.powerStart;
			end
			x_shifted = x - powerStart_h;
			if x_shifted<0
				s = 1;
			else
				s = obj.positiveScaling;
			end
		end
		
	
		
		
		
	end

end
