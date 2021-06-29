%
%  Template for experiment files. Provide here general information on your
%  simulation study, e.g.: This file contains experiments for
%  spatio-temporal prediction of time series with drilling data. The
%  experiments 1002, 3004, 4003 and 4004 correspond to the paper ...
%
%

classdef UAVPlacementExperiments < ExperimentFunctionSet
	
	properties
		% You may define here constants that you will use in the
		% experiments of this file
		
	end
	
	methods
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 10. Experiments in 1D
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% Investigate the choice of the utility function in terms of the optimal
		% spatial distribution of the BSs.
		% - Infinitely many MUs
		% Observations: there is an optimum distance. 
		function F = experiment_1001(obj,niter)
			
			v_txPowerBSs = [2;3];
			v_limitsRegion = [0,10];
			positionBS1 = 5;			
				
			gainConstant = 1;
			height = 1;
			f_gain = @(xBS,xMU) gainConstant./((xBS-xMU).^2+height^2); % propagation gain for distance d
			f_gainDXBS = @(xBS,xMU) 2*f_gain(xBS,xMU).^2/gainConstant.*(xMU-xBS); % derivative w.r.t. xBS
			
			utilityFunctionSelector = 3;
			if utilityFunctionSelector == 1 % utility = capacity
				ch_utilityName = 'Capacity';
				noisePower = 1;
				f_utility = @(p) log2(1+p/noisePower);
				f_utilityPrime = @(p) 1./(1+p/noisePower)/noisePower/log(2);
			elseif utilityFunctionSelector == 2 % sigmoid
				ch_utilityName = 'Sigmoid';
				min_rxPower = 1;
				max_rxPower = 2;
				f_sigmoid = @(x) exp(x)./(exp(x)+1);
				f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
				f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
				f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
					* 6/(max_rxPower-min_rxPower);
			elseif utilityFunctionSelector == 3 % leaky ReLU
				ch_utilityName = 'Leaky ReLU';
				powerStart = 1;
				positiveScaling = 0.0002;
				f_utility = @(p) min( p - powerStart ,positiveScaling*(p - powerStart));
				f_utilityPrime = @(p) positiveScaling*(p>powerStart) + 1*(p<=powerStart); 
			end
						
			v_positionsMUs = linspace(v_limitsRegion(1),v_limitsRegion(2),200)'; % for displaying purposes and estimating expectations
			
			% Figure with power and utility averages for a small set of BS2
			% locations.
			v_positionsBS2 = linspace(v_limitsRegion(1),mean(v_limitsRegion),6);
			v_powerBS1 = v_txPowerBSs(1)*f_gain(positionBS1,v_positionsMUs);
			for ind_positionBS2 = 1:length(v_positionsBS2)								
				v_powerBS2 = v_txPowerBSs(2)*f_gain(v_positionsBS2(ind_positionBS2),v_positionsMUs);
				v_utility = f_utility( v_powerBS1 + v_powerBS2 );
				ch_title = sprintf('Avg. utility = %g',mean(v_utility));				
				m_multiplot(1,ind_positionBS2) = GFigure('m_X',v_positionsMUs','m_Y',GFigure.formMatrixByPaddingRows(...
					v_powerBS1'+v_powerBS2',v_powerBS1',v_powerBS2'),'c_legend',{'Total power','Power BS1',...
					'Power BS2'},'ch_xlabel','Location','ch_ylabel','Power','ch_title',ch_title);				
				m_multiplot(2,ind_positionBS2) = GFigure('m_X',v_positionsMUs','m_Y',v_utility',...
					'ch_xlabel','Location','ch_ylabel','Utility');				
			end
			F(1) = GFigure('m_multiplot',m_multiplot);
			
			% Figure with average utility and average gradient
			num_points = 20;
			m_multiplot = GFigure();
			v_positionsBS2 = linspace(v_limitsRegion(1),mean(v_limitsRegion),num_points);
			v_powerBS1 = v_txPowerBSs(1)*f_gain(positionBS1,v_positionsMUs);		
			for ind_positionBS2 = length(v_positionsBS2):-1:1
				v_powerBS2 = v_txPowerBSs(2)*f_gain(v_positionsBS2(ind_positionBS2),v_positionsMUs);
				v_utility = f_utility( v_powerBS1 + v_powerBS2 );
				v_avgUtility(ind_positionBS2) = mean(v_utility);
				
				v_gradient = f_utilityPrime( v_powerBS1 + v_powerBS2 )...
					.* v_txPowerBSs(2).* f_gainDXBS( v_positionsBS2(ind_positionBS2) , v_positionsMUs );
				v_avgGradient(ind_positionBS2) = mean(v_gradient);						
			end
			m_multiplot(1,ind_positionBS2) = GFigure('m_X',v_positionsBS2,'m_Y',v_avgUtility,...
				'ch_xlabel','Location','ch_ylabel','Average Utility',...
				'ch_title',sprintf('Utility = %s',ch_utilityName));
			m_multiplot(2,ind_positionBS2) = GFigure('m_X',v_positionsBS2,'m_Y',v_avgGradient,...
				'ch_xlabel','Location','ch_ylabel','Average Gradient');							
			F(2) = GFigure('m_multiplot',m_multiplot);
			
		end	
	
		
		% Aerial BSs seek their optimum location using gradient descent
		% - Infinitely many MUs
		% Observation: They indeed spread nicely across the space. 
		function F = experiment_1002(obj,niter)
				
			v_txPowerBSs = [2;3;4];
			v_limitsRegion = [0,10];
			v_initialPositionsBS = [4;5;10];
			v_stepSizes = [0.5;0.5;0.5]; % step size for each BS
				
			gainConstant = 1;
			height = 1;
			f_gain = @(xBS,xMU) gainConstant./((xBS-xMU).^2+height^2); % propagation gain for distance d
			f_gainDXBS = @(xBS,xMU) 2*f_gain(xBS,xMU).^2/gainConstant.*(xMU-xBS); % derivative w.r.t. xBS
						
			min_rxPower = 1;
			max_rxPower = 2;
			f_sigmoid = @(x) exp(x)./(exp(x)+1);
			f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
			f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
			f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
				* 6/(max_rxPower-min_rxPower);
									
			v_positionsMUs = linspace(v_limitsRegion(1),v_limitsRegion(2),200); % for displaying purposes and estimating expectations
			
			% Figure with power and utility function
			v_positionsBSs = v_initialPositionsBS;
			num_BS = length(v_positionsBSs);
			num_MU = length(v_positionsMUs);	
			for ind_BS = num_BS+1:-1:2
				c_legend{ind_BS} = sprintf('Power BS%d',ind_BS);
			end
			c_legend{1} = 'Total power';
			while 1
				
				% Power and utility
				m_power = diag(v_txPowerBSs)*f_gain(repmat(v_positionsBSs,1,num_MU),...
					repmat(v_positionsMUs,num_BS,1));	
				v_totalPower1D = sum(m_power,1);
				v_utility = f_utility( v_totalPower1D );					
				
				% Gradient
				m_gradient = repmat(f_utilityPrime( v_totalPower1D ),num_BS,1)...
					.* (diag(v_txPowerBSs)* f_gainDXBS( repmat(v_positionsBSs,1,num_MU) , ...
					repmat(v_positionsMUs,num_BS,1)));%      one row per BS, one col per MU
				v_avgGradient = mean(m_gradient,2);
												
				ch_title = sprintf('Avg. utility = %g, Avg. grad = %s',mean(v_utility),vector_to_str(v_avgGradient));				
				m_multiplot(1,1) = GFigure('m_X',v_positionsMUs,'m_Y',[v_totalPower1D;m_power],...
					'c_legend',c_legend,'ch_xlabel','Location','ch_ylabel','Power','ch_title',ch_title);				
				m_multiplot(2,1) = GFigure('m_X',v_positionsMUs,'m_Y',v_utility,...
					'ch_xlabel','Location','ch_ylabel','Utility');				
				F(1) = GFigure('m_multiplot',m_multiplot);
				F.plot();
				pause()
				
				% Gradient step
				v_positionsBSs = v_positionsBSs + diag(v_stepSizes)*v_avgGradient;
				
			end
			
			
		end	
	
		
		% Aerial BSs seek their optimum location using STOCHASTIC gradient
		% descent with minibatches.
		% - Infinitely many MUs in the sense that the MUs reporting
		% the control information at every time instant are chosen
		% according to a uniform distribution among a continuous set of MUs
		% on the real line. 
		% Observation: They indeed spread nicely across the space. 
		function F = experiment_1003(obj,niter)
				
			v_txPowerBSs = [2;3;4;6];
			v_limitsRegion = [0,10];
			v_initialPositionsBS = [10;10.1;10.2;10.3]; %[4;5;10;10.1];
			v_stepSizes = 1*ones(length(v_txPowerBSs),1); % step size for each BS
			miniBatchSize = 5; % number of measurements per update
				
			gainConstant = 1;
			height = 1;
			f_gain = @(xBS,xMU) gainConstant./((xBS-xMU).^2+height^2); % propagation gain for distance d
			f_gainDXBS = @(xBS,xMU) 2*f_gain(xBS,xMU).^2/gainConstant.*(xMU-xBS); % derivative w.r.t. xBS
						
% 			min_rxPower = 2;
% 			max_rxPower = 3;
% 			f_sigmoid = @(x) exp(x)./(exp(x)+1);
% 			f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
% 			f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
% 			f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
% 				* 6/(max_rxPower-min_rxPower);
			
			utilityFunctionSelector = 3;
			if utilityFunctionSelector == 1 % utility = capacity
				ch_utilityName = 'Capacity';
				noisePower = 1;
				f_utility = @(p) log2(1+p/noisePower);
				f_utilityPrime = @(p) 1./(1+p/noisePower)/noisePower/log(2);
			elseif utilityFunctionSelector == 2 % sigmoid
				ch_utilityName = 'Sigmoid';
				min_rxPower = 1;
				max_rxPower = 2;
				f_sigmoid = @(x) exp(x)./(exp(x)+1);
				f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
				f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
				f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
					* 6/(max_rxPower-min_rxPower);
			elseif utilityFunctionSelector == 3 % leaky ReLU
				ch_utilityName = 'Leaky ReLU';
				powerStart = 2;
				positiveScaling = 0.0002;
				f_utility = @(p) min( p - powerStart ,positiveScaling*(p - powerStart));
				f_utilityPrime = @(p) positiveScaling*(p>powerStart) + 1*(p<=powerStart); 
			end
			
			v_positionsMUTest = linspace(v_limitsRegion(1),v_limitsRegion(2),100); % for displaying purposes 
			
			% Figure with power and utility function
			v_positionsBSs = v_initialPositionsBS;
			num_BS = length(v_positionsBSs);
			for ind_BS = num_BS+2:-1:3
				c_legend{ind_BS} = sprintf('Power BS%d',ind_BS);
				c_styles{ind_BS} = '--';
			end
			c_legend{1} = 'MU sending measurements';
			c_styles{1} = 'x';
			c_legend{2} = 'Total power';			
			c_styles{2} = '-';
			while 1
				
				% Power and utility
				[v_totalPower1D , m_partialPower] = UAVPlacementExperiments.totalPower1D( ...
					v_txPowerBSs , v_positionsBSs , f_gain , v_positionsMUTest );				
				v_utilityTest = f_utility( v_totalPower1D );					
				
				% Gradient
				v_positionsMUMeasurement = v_limitsRegion(1) + (v_limitsRegion(2)-v_limitsRegion(1))*...
					rand(1,miniBatchSize);
				[v_avgGradient] = UAVPlacementExperiments.stochasticGradient1D( ...
					v_txPowerBSs , v_positionsBSs , f_gain , f_gainDXBS, f_utilityPrime , v_positionsMUMeasurement );
												
				ch_title = sprintf('Avg. utility = %g, Avg. grad = %s',mean(v_utilityTest),vector_to_str(v_avgGradient));				
				m_multiplot(1,1) = GFigure('m_X',...
					GFigure.formMatrixByPaddingRows(v_positionsMUMeasurement,repmat(v_positionsMUTest,num_BS+1,1)),...
					'm_Y',GFigure.formMatrixByPaddingRows(zeros(1,miniBatchSize),v_totalPower1D,m_partialPower),...
					'c_legend',c_legend,'ch_xlabel','Location','ch_ylabel','Power','ch_title',ch_title,...
					'c_styles',c_styles);				
				m_multiplot(2,1) = GFigure('m_X',v_positionsMUTest,'m_Y',v_utilityTest,...
					'ch_xlabel','Location','ch_ylabel','Utility');				
				F(1) = GFigure('m_multiplot',m_multiplot);
				F.plot();
				pause(.05)
				
				% Gradient step
				v_positionsBSs = v_positionsBSs + diag(v_stepSizes)*v_avgGradient;
				
			end
			
			
		end	
	
		
		% Aerial BSs seek their optimum location using stochastic gradient
		% descent with minibatches.
		% - FINITELY many MUs. All same throughput
		function F = experiment_1004(obj,niter)
				
			% Base stations
			v_txPowerBSs = [2;3;4;6];
			v_limitsRegion = [0,10];
			v_initialPositionsBS = [4;5;10;10.1];
			v_stepSizes = 1*ones(length(v_txPowerBSs),1); % step size for each BS
			height = 1;
			
			% Propagation
			gainConstant = 1;			
			f_gain = @(xBS,xMU) gainConstant./((xBS-xMU).^2+height^2); % propagation gain for distance d
			f_gainDXBS = @(xBS,xMU) 2*f_gain(xBS,xMU).^2/gainConstant.*(xMU-xBS); % derivative w.r.t. xBS
						
			% Mobile users
			miniBatchSize = 5; % number of measurements per update
			num_MU = 40;
			v_positionsMU = v_limitsRegion(1) + (v_limitsRegion(2)-v_limitsRegion(1))*...
					rand(1,num_MU); % where the utility is measured
			
			% Target received levels			
			utilityFunctionSelector = 3;
			if utilityFunctionSelector == 1 % utility = capacity
				ch_utilityName = 'Capacity';
				noisePower = 1;
				f_utility = @(p) log2(1+p/noisePower);
				f_utilityPrime = @(p) 1./(1+p/noisePower)/noisePower/log(2);
			elseif utilityFunctionSelector == 2 % sigmoid
				ch_utilityName = 'Sigmoid';
				min_rxPower = 2;
				max_rxPower = 3;
				f_sigmoid = @(x) exp(x)./(exp(x)+1);
				f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
				f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
				f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
					* 6/(max_rxPower-min_rxPower);
			elseif utilityFunctionSelector == 3 % leaky ReLU
				ch_utilityName = 'Leaky ReLU';
				powerStart = 2;
				positiveScaling = 0.0002;
				f_utility = @(p) min( p - powerStart ,positiveScaling*(p - powerStart));
				f_utilityPrime = @(p) positiveScaling*(p>powerStart) + 1*(p<=powerStart); 
			end
			
			% Figure with power and utility function
			v_positionsBSs = v_initialPositionsBS;
			v_positionsMUPlot = linspace(v_limitsRegion(1),v_limitsRegion(2),100); % where the power is plotted
			num_BS = length(v_positionsBSs);
			for ind_BS = num_BS+3:-1:4
				c_legend{ind_BS} = sprintf('Power BS%d',ind_BS);
				c_styles{ind_BS} = '--';
			end
			c_legend{1} = 'MUs';
			c_styles{1} = 'o';
			c_legend{2} = 'MUs sending measurements';
			c_styles{2} = 'x';
			c_legend{3} = 'Total power';			
			c_styles{3} = '-';
			while 1
				
				% Power and utility for displaying purposes
				[v_totalPower1D , m_partialPower] = UAVPlacementExperiments.totalPower1D( ...
					v_txPowerBSs , v_positionsBSs , f_gain , v_positionsMUPlot );				
				v_totalPower1DMUs = UAVPlacementExperiments.totalPower1D( ...
					v_txPowerBSs , v_positionsBSs , f_gain , v_positionsMU );							
				v_utilityTest = f_utility( v_totalPower1DMUs );					
				
				% Gradient
				v_positionsMUMeasurement = v_positionsMU( randperm( num_MU , miniBatchSize ) );				
				[v_avgGradient] = UAVPlacementExperiments.stochasticGradient1D( ...
					v_txPowerBSs , v_positionsBSs , f_gain , f_gainDXBS, f_utilityPrime , v_positionsMUMeasurement );
												
				ch_title = sprintf('Avg. utility = %g, Avg. grad = %s',mean(v_utilityTest),vector_to_str(v_avgGradient));				
				m_multiplot(1,1) = GFigure('m_X',...
					GFigure.formMatrixByPaddingRows(v_positionsMU,v_positionsMUMeasurement,repmat(v_positionsMUPlot,num_BS+1,1)),...
					'm_Y',GFigure.formMatrixByPaddingRows(zeros(1,num_MU),zeros(1,miniBatchSize),v_totalPower1D,m_partialPower),...
					'c_legend',c_legend,'ch_xlabel','Location','ch_ylabel','Power','ch_title',ch_title,...
					'c_styles',c_styles);				
				m_multiplot(2,1) = GFigure('m_X',v_positionsMU,'m_Y',v_utilityTest,...
					'ch_xlabel','Location','ch_ylabel','Utility','ch_plotType2D','stem');				
				F(1) = GFigure('m_multiplot',m_multiplot);
				F.plot();
				pause(.05)
				
				% Gradient step
				v_positionsBSs = v_positionsBSs + diag(v_stepSizes)*v_avgGradient;
				
			end
			
			
		end	
	
		% Gradient and utility for the position of an AirBS when the other
		% AirBSs are fixed.
		% TO DO: 
		%  - in a Different figure -> 
		%       -Plot power of each non-test base station across space
		%       -plot the utility for the MUs when the test
		%        BS is not present
		function F = experiment_1005(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			displaySettings.v_windowPosition = [   744-500        449        1094         601];
			
			rng(1)
			
			regionSideLength = 7; % km % region is a line segment
			wattPerPowerUnit = 1; %3.971641173621411e-13; % scaling to 
			%wattPerPowerUnit = 3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 2;
			b_2D = 0;
			%powerStart = 2;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','leakyRelu','ch_aggregate','smoothMaxPower',...
			 	'powerStart',powerStart , 'b_dBm',1 );
			%templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','totalPower',...
			%	'powerStart',powerStart ,'powerStop',powerStop, 'b_dBm',1, 'positiveScaling',0.0 );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1,b_2D);
				
			% Base stations
			height = 0.03;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',1e4,... %25
				'v_position',[0,0,height]','b_fixedHeight',1,'b_debug',1);
			%v_txPowerBSs = dBm_to_Watt([7,9,9,9,10,10,12]-3)/wattPerPowerUnit;
			v_txPowerBSs = dBm_to_Watt([6,4])/wattPerPowerUnit;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength,1,b_2D);
			
			% Channel
			f=2390e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);
			channel = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plot of utility and gradient as a function of the position of
			% v_aerialBaseStations(1)			
			v_txPower = v_aerialBaseStations.powerAerialBaseStations();
			v_xaxisTestPositions = linspace(-regionSideLength,2*regionSideLength,1000);
			v_gradientWrtXPosition = zeros(1,length(v_xaxisTestPositions));
			for ind_point = length(v_xaxisTestPositions):-1:1
				v_aerialBaseStations(1).v_position(1) = v_xaxisTestPositions(ind_point);
				m_positionsBSs = v_aerialBaseStations.positionAerialBaseStations();
				
				for ind_selectedMU = length(v_mobileUsers):-1:1
					
					selectedMU = v_mobileUsers(ind_selectedMU);
					% Utility gradients w.r.t. power from MUs
					v_rxPower = v_txPower.*channel.gain(m_positionsBSs,selectedMU.v_position);
					v_utilityGradientWrtRxPowerByMU(:,ind_selectedMU) = selectedMU.utilityGradient(v_rxPower);
					v_utilityByMU(1,ind_selectedMU) = selectedMU.utility(v_rxPower);
					%    One row per base station, one column per MU
					
					v_RxPowerGradientWrtXPosition = v_aerialBaseStations(1).power * channel.gainGradient(v_aerialBaseStations(1).v_position,selectedMU.v_position);
					
					v_gradientWrtXPosition(ind_point) = v_gradientWrtXPosition(ind_point) + v_RxPowerGradientWrtXPosition(1)*v_utilityGradientWrtRxPowerByMU(ind_selectedMU);
					if ind_selectedMU == 1
						v_rxPowerMU1FromBS1(ind_point) = v_aerialBaseStations(1).power * channel.gain(v_aerialBaseStations(1).v_position,selectedMU.v_position);
						v_rxPowerGradientMU1FromBS1(:,ind_point) = v_aerialBaseStations(1).power * channel.gainGradient(v_aerialBaseStations(1).v_position,selectedMU.v_position);
					end
				end
				v_utility(1,ind_point) = mean(v_utilityByMU,2);
				v_utilityGradientWrtRxPowerBS1(:,ind_point) = sum(v_utilityGradientWrtRxPowerByMU(1,:),2);		
			end
			
			
			m_multiplot(1,1) = GFigure('m_X',v_xaxisTestPositions,'m_Y',v_utility,...
				'ch_xlabel','X-Position of BS1 [km]','ch_ylabel','Utility','ch_title','Utility for the position of BS1');
			m_positionsMUs = v_mobileUsers.positionMobileUsers();
			m_multiplot(1,1,2) = GFigure('m_X',m_positionsMUs(1,:),'m_Y',mean(v_utility)*ones(size(m_positionsMUs(1,:))),...
				'c_styles',{'x'},'m_colorset',[1 0 0]);
			m_multiplot(1,1,3) = GFigure('m_X',m_positionsBSs(1,2:end),'m_Y',mean(v_utility)*ones(1,length(v_txPowerBSs)-1),...
				'c_styles',{'o'});
			m_multiplot(1,1).c_legend = {'curve','MU positions','BS positions'};
			
			m_multiplot(2,1) = GFigure('m_X',v_xaxisTestPositions,'m_Y',v_utilityGradientWrtRxPowerBS1,...
				'ch_xlabel','X-Position of BS1 [km]','ch_ylabel','Gradient','ch_title','Utility Gradient Wrt Power Rx from BS1');
			m_multiplot(3,1) = GFigure('m_X',v_xaxisTestPositions,'m_Y',v_gradientWrtXPosition,...
				'ch_xlabel','X-Position of BS1 [km]','ch_ylabel','Gradient','ch_title','Utility Gradient wrt position of BS1');		
			F(1) = GFigure('m_multiplot',m_multiplot);
			
			m_multiplot(1,1) = GFigure('m_X',v_xaxisTestPositions,'m_Y',v_rxPowerMU1FromBS1,...
				'ch_xlabel','X-Position of BS1 [km]','ch_ylabel','Utility','ch_title','Power from BS1 rx at MU1');
			m_multiplot(2,1) = GFigure('m_X',v_xaxisTestPositions,'m_Y',v_rxPowerGradientMU1FromBS1(1,:),...
				'ch_xlabel','X-Position of BS1 [km]','ch_ylabel','Utility','ch_title','Gradient of Power from BS1 rx at MU1 wrt position of BS1');
			%F(2) = GFigure('m_multiplot',m_multiplot(1:2,1));
						
		end	
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 20. Experiments in 2D
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% Aerial BSs seek their optimum location using stochastic gradient
		% descent with minibatches.
		% - finitely many MUs. All same throughput.
		% Observations: the BSs nicely spread across space.
		function F = experiment_2001(obj,niter)
				
			global displaySettings 
			displaySettings.v_windowPosition = [   744         449        1094         601];
			num_iter = 100;
			
			% Base stations
			v_txPowerBSs = [2,3,3,3,4,4,6];
			num_BS = length(v_txPowerBSs);
			regionSideLength = 7; % region is a square with bottom left edge at the origin. This is the side
			m_initialPositionsBS = regionSideLength/3*rand(2,num_BS);%[ 2,2;3,8;2,8;9,1]'; % one column per BS
			v_stepSizes = 25*ones(length(v_txPowerBSs),1); % step size for each BS
			height = 1;
			
			% Propagation
			gainConstant = 1;			
			f_gain = @(v_xBS,v_xMU) gainConstant./(sum((v_xBS-v_xMU).^2,1)+height^2); % propagation gain for distance d
			f_gainDXBS = @(v_xBS,v_xMU) 2/gainConstant*(v_xMU-v_xBS).*f_gain(v_xBS,v_xMU).^2; % gradient w.r.t. v_xBS
						
			% Mobile users
			miniBatchSize = 200; % number of measurements per update
			num_MU = 200;
			m_positionsMU = regionSideLength*rand(2,num_MU); % where the utility is measured
			
			% Target received levels			
			utilityFunctionSelector = 3;
			if utilityFunctionSelector == 1 % utility = capacity
				ch_utilityName = 'Capacity';
				noisePower = 1;
				f_utility = @(p) log2(1+p/noisePower);
				f_utilityPrime = @(p) 1./(1+p/noisePower)/noisePower/log(2);
			elseif utilityFunctionSelector == 2 % sigmoid
				ch_utilityName = 'Sigmoid';
				min_rxPower = 1;
				max_rxPower = 2;
				f_sigmoid = @(x) exp(x)./(exp(x)+1);
				f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
				f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
				f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
					* 6/(max_rxPower-min_rxPower);
			elseif utilityFunctionSelector == 3 % leaky ReLU
				ch_utilityName = 'Leaky ReLU';
				powerStart = 3;
				min_rxPower = powerStart;
				max_rxPower = powerStart+1;
				positiveScaling = 0.0002;
				f_utility = @(p) min( p - powerStart ,positiveScaling*(p - powerStart));
				f_utilityPrime = @(p) positiveScaling*(p>powerStart) + 1*(p<=powerStart); 
			end
			
			% Figure with power and utility function
			m_positionsBSs = m_initialPositionsBS;
			num_pointsGridPerDim = 30;
			[m_xCoordinatesMUPlot,m_yCoordinatesMUPlot] = meshgrid(...
				linspace(0,regionSideLength,num_pointsGridPerDim),...
				linspace(0,regionSideLength,num_pointsGridPerDim)); % where the power is plotted
			%m_yCoordinatesMUPlot = flipud(m_yCoordinatesMUPlot); % Y axis increases upwards
			m_positionsMUPlot = [m_xCoordinatesMUPlot(:),m_yCoordinatesMUPlot(:)]';
			v_fractionUsersOverMin_rxPower = [];
			v_fractionUsersOverMax_rxPower = [];
			v_fractionUsersWithinPowerInterval = [];
			
			for ind_iter = 1:num_iter
				
				% Power and utility for displaying purposes
				[v_totalPower ] = UAVPlacementExperiments.totalPower2D( ...
					v_txPowerBSs , m_positionsBSs , f_gain , m_positionsMUPlot );	
				v_utilityEverywhere = f_utility( v_totalPower );
				v_totalPowerMUs = UAVPlacementExperiments.totalPower2D( ...
					v_txPowerBSs , m_positionsBSs , f_gain , m_positionsMU );							
				v_utilityMU = f_utility( v_totalPowerMUs );					
						
				% Fraction of users with power > min_rxPower vs time 
				v_fractionUsersOverMin_rxPower(end+1) = mean( v_totalPowerMUs > min_rxPower );
				v_fractionUsersOverMax_rxPower(end+1) = mean( v_totalPowerMUs > max_rxPower );
				v_fractionUsersWithinPowerInterval(end+1) = mean( (v_totalPowerMUs <= max_rxPower )...
					& (v_totalPowerMUs >= min_rxPower ) );
				
	            % Histogram of SNR for unit noise power
				v_binCenters = linspace(0,max(v_totalPower),15);
				[v_num_MUEachBin] = hist( v_totalPowerMUs , v_binCenters);
				v_pdfPower = v_num_MUEachBin/num_MU/(v_binCenters(2)-v_binCenters(1));
				
				% Figures
				if ind_iter == 1
					%v_caxisPower = [0,.5*max(v_totalPower)];
					v_caxisPower = [min_rxPower,max_rxPower];
				end
				ch_title = sprintf('Avg. utility = %g',mean(v_utilityMU));
				m_multiplot(1,1,1) = GFigure('m_X',m_positionsMU(1,:),'m_Y',m_positionsMU(2,:),...
					'm_Z',zeros(1,num_MU),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Height',...
					'ch_plotType3D','stem3','c_styles',{'.k'},...
					'ch_title','Locations of MUs and Aerial BSs');
					%'c_legend',{'Mobile Users','Aerial BSs'});		
				m_multiplot(1,1,2) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
					'm_Z',height*ones(1,num_BS),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_plotType3D','stem3','c_styles',{'sb'});		
				m_multiplot(1,2) = GFigure('m_X',m_xCoordinatesMUPlot,...
					'm_Y',m_yCoordinatesMUPlot,...
					'm_Z',reshape(v_totalPower,num_pointsGridPerDim,num_pointsGridPerDim),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Total power',...					
					'ch_plotType3D','surf');			
% 				m_multiplot(1,2) = GFigure('m_X',m_xCoordinatesMUPlot,...
% 					'm_Y',m_yCoordinatesMUPlot,...
% 					'm_Z',reshape(v_utilityEverywhere,num_pointsGridPerDim,num_pointsGridPerDim),...
% 					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...					
% 					'ch_zlabel','Utility function',...
% 					'ch_plotType3D','surf');
				m_multiplot(1,3,1) = GFigure('m_X',1:length(v_fractionUsersOverMin_rxPower),...
					'm_Y',[v_fractionUsersOverMin_rxPower;v_fractionUsersOverMax_rxPower;...
					v_fractionUsersWithinPowerInterval],...
					'ch_xlabel','Iteration index',...
					'c_legend',{'Frac. MUs over Pmin','Frac. MUs over Pmax',...
					'Frac. MUs within interval'},'ch_legendPosition','northwest');
				m_multiplot(2,1,1) = GFigure('m_X',m_xCoordinatesMUPlot,...
					'm_Y',m_yCoordinatesMUPlot,...
					'm_Z',reshape(v_totalPower,num_pointsGridPerDim,num_pointsGridPerDim),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Total power',...
					'ch_plotType3D','imagesc','v_caxis',v_caxisPower,...
					'ch_title','Power');	
				m_multiplot(2,1,2) = GFigure('m_X',m_positionsMU(1,:),'m_Y',m_positionsMU(2,:),...
					'm_Z',zeros(1,num_MU),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_plotType3D','stem3','c_styles',{'.w'});		
				m_multiplot(2,1,3) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
					'm_Z',zeros(1,num_BS),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_plotType3D','stem3','c_styles',{'sc'});						
				m_multiplot(2,2) = GFigure('m_X',m_positionsMU(1,:),'m_Y',m_positionsMU(2,:),...
					'm_Z',v_utilityMU,...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Utility',...
					'ch_title',ch_title,...
					'ch_zlabel','Utility','ch_plotType3D','stem3'); %,'c_styles','filled');				
				m_multiplot(2,3,1) = GFigure('m_X',v_binCenters,...
					'm_Y',v_pdfPower,...
					'c_legend',{'Distribution','Pmin','Pmax'},...
					'ch_xlabel','Received power','v_ylim',[0,1.2*max(v_pdfPower)],...
					'ch_ylabel','Distribution MU Rx Power');
				m_multiplot(2,3,2) = GFigure('m_X',	[min_rxPower,min_rxPower;max_rxPower,max_rxPower],...
					'm_Y',[0,1.2*max(v_pdfPower);0,1.2*max(v_pdfPower)]);									
				F(1) = GFigure('m_multiplot',m_multiplot);
				F.plot();
				v_frames(ind_iter) = getframe(gcf);
				%pause();
				pause(.005)
				
				
				% Gradient
				v_positionsMUMeasurement = m_positionsMU( : , randperm( num_MU , miniBatchSize ) );				
				[m_avgGradient] = UAVPlacementExperiments.stochasticGradient2D( ...
					v_txPowerBSs , m_positionsBSs , f_gain , f_gainDXBS, f_utilityPrime , v_positionsMUMeasurement );				
								
				% Gradient step
				m_positionsBSs = m_positionsBSs + m_avgGradient*diag(v_stepSizes);
				
			end
			
			
			% Create video
			writerObj = VideoWriter('./Experiments/UAVPlacementExperiments_data/aerialBaseStations',...
				'Motion JPEG AVI');
		%	'Motion JPEG 2000');
			writerObj.FrameRate = 10;
			writerObj.open();
			for ind_frame = 1:length(v_frames)
				writerObj.writeVideo(v_frames(ind_frame));				
			end
			writerObj.close();
			
		end	
		
		% same as 2001 but without showing Pmax
		function F = experiment_2002(obj,niter)
				
			global displaySettings 
			displaySettings.v_windowPosition = [   744         449        1094         601];
			num_iter = 100;
			
			% Base stations
			v_txPowerBSs = [2,3,3,3,4,4,6];
			num_BS = length(v_txPowerBSs);
			regionSideLength = 7; % region is a square with bottom left edge at the origin. This is the side
			m_initialPositionsBS = regionSideLength/3*rand(2,num_BS);%[ 2,2;3,8;2,8;9,1]'; % one column per BS
			v_stepSizes = 15*ones(length(v_txPowerBSs),1); % step size for each BS
			height = 1;
			
			% Propagation
			gainConstant = 1;			
			f_gain = @(v_xBS,v_xMU) gainConstant./(sum((v_xBS-v_xMU).^2,1)+height^2); % propagation gain for distance d
			f_gainDXBS = @(v_xBS,v_xMU) 2/gainConstant*(v_xMU-v_xBS).*f_gain(v_xBS,v_xMU).^2; % gradient w.r.t. v_xBS
						
			% Mobile users
			miniBatchSize = 200; % number of measurements per update
			num_MU = 200;
			m_positionsMU = regionSideLength*rand(2,num_MU); % where the utility is measured
			
			% Target received levels			
			min_rxPower = 2;
			max_rxPower = 3;
			f_sigmoid = @(x) exp(x)./(exp(x)+1);
			f_sigmoidPrime = @(x) exp(-x)./(1+exp(-x)).^2;
			f_utility = @(p) f_sigmoid(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower));
			f_utilityPrime = @(p) f_sigmoidPrime(3*(2*p-min_rxPower-max_rxPower)/(max_rxPower-min_rxPower))...
				* 6/(max_rxPower-min_rxPower);
			
			% Figure with power and utility function
			m_positionsBSs = m_initialPositionsBS;
			num_pointsGridPerDim = 30;
			[m_xCoordinatesMUPlot,m_yCoordinatesMUPlot] = meshgrid(...
				linspace(0,regionSideLength,num_pointsGridPerDim),...
				linspace(0,regionSideLength,num_pointsGridPerDim)); % where the power is plotted
			%m_yCoordinatesMUPlot = flipud(m_yCoordinatesMUPlot); % Y axis increases upwards
			m_positionsMUPlot = [m_xCoordinatesMUPlot(:),m_yCoordinatesMUPlot(:)]';
			v_fractionUsersOverMin_rxPower = [];
			v_fractionUsersOverMax_rxPower = [];
			v_fractionUsersWithinPowerInterval = [];
			
			for ind_iter = 1:num_iter
				
				% Power and utility for displaying purposes
				[v_totalPower ] = UAVPlacementExperiments.totalPower2D( ...
					v_txPowerBSs , m_positionsBSs , f_gain , m_positionsMUPlot );	
				v_utilityEverywhere = f_utility( v_totalPower );
				v_totalPowerMUs = UAVPlacementExperiments.totalPower2D( ...
					v_txPowerBSs , m_positionsBSs , f_gain , m_positionsMU );							
				v_utilityMU = f_utility( v_totalPowerMUs );					
						
				% Fraction of users with power > min_rxPower vs time 
				v_fractionUsersOverMin_rxPower(end+1) = mean( v_totalPowerMUs > min_rxPower );
				v_fractionUsersOverMax_rxPower(end+1) = mean( v_totalPowerMUs > max_rxPower );
				v_fractionUsersWithinPowerInterval(end+1) = mean( (v_totalPowerMUs <= max_rxPower )...
					& (v_totalPowerMUs >= min_rxPower ) );
				
	            % Histogram of SNR for unit noise power
				v_binCenters = linspace(0,max(v_totalPower),15);
				[v_num_MUEachBin] = hist( v_totalPowerMUs , v_binCenters);
				v_pdfPower = v_num_MUEachBin/num_MU/(v_binCenters(2)-v_binCenters(1));
				
				% Figures
				if ind_iter == 1
					%v_caxisPower = [0,.5*max(v_totalPower)];
					v_caxisPower = [min_rxPower,max_rxPower];
				end
				ch_title = sprintf('Avg. utility = %g',mean(v_utilityMU));
				m_multiplot(1,1,1) = GFigure('m_X',m_positionsMU(1,:),'m_Y',m_positionsMU(2,:),...
					'm_Z',zeros(1,num_MU),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Height',...
					'ch_plotType3D','stem3','c_styles',{'.k'},...
					'ch_title','Locations of MUs and Aerial BSs');
					%'c_legend',{'Mobile Users','Aerial BSs'});		
				m_multiplot(1,1,2) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
					'm_Z',height*ones(1,num_BS),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_plotType3D','stem3','c_styles',{'sb'});		
				m_multiplot(1,2) = GFigure('m_X',m_xCoordinatesMUPlot,...
					'm_Y',m_yCoordinatesMUPlot,...
					'm_Z',reshape(v_totalPower,num_pointsGridPerDim,num_pointsGridPerDim),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Total power',...					
					'ch_plotType3D','surf');			
% 				m_multiplot(1,2) = GFigure('m_X',m_xCoordinatesMUPlot,...
% 					'm_Y',m_yCoordinatesMUPlot,...
% 					'm_Z',reshape(v_utilityEverywhere,num_pointsGridPerDim,num_pointsGridPerDim),...
% 					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...					
% 					'ch_zlabel','Utility function',...
% 					'ch_plotType3D','surf');
				m_multiplot(1,3,1) = GFigure('m_X',1:length(v_fractionUsersOverMin_rxPower),...
					'm_Y',[v_fractionUsersOverMin_rxPower],...
					'ch_xlabel','Iteration index',...
					'c_legend',{'Frac. MUs over Pmin'},'ch_legendPosition','northwest',...
					'v_ylim',[0,1]);
				m_multiplot(2,1,1) = GFigure('m_X',m_xCoordinatesMUPlot,...
					'm_Y',m_yCoordinatesMUPlot,...
					'm_Z',reshape(v_totalPower,num_pointsGridPerDim,num_pointsGridPerDim),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Total power',...
					'ch_plotType3D','imagesc','v_caxis',v_caxisPower,...
					'ch_title','Power');	
				m_multiplot(2,1,2) = GFigure('m_X',m_positionsMU(1,:),'m_Y',m_positionsMU(2,:),...
					'm_Z',zeros(1,num_MU),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_plotType3D','stem3','c_styles',{'.w'});		
				m_multiplot(2,1,3) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
					'm_Z',zeros(1,num_BS),...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_plotType3D','stem3','c_styles',{'sc'});						
				m_multiplot(2,2) = GFigure('m_X',m_positionsMU(1,:),'m_Y',m_positionsMU(2,:),...
					'm_Z',v_utilityMU,...
					'ch_xlabel','x coordinate','ch_ylabel','y coordinate',...
					'ch_zlabel','Utility',...
					'ch_title',ch_title,...
					'ch_zlabel','Utility','ch_plotType3D','stem3'); %,'c_styles','filled');				
				m_multiplot(2,3,1) = GFigure('m_X',v_binCenters,...
					'm_Y',v_pdfPower,...
					'c_legend',{'Distribution','Pmin'},...
					'ch_xlabel','Received power','v_ylim',[0,1.2*max(v_pdfPower)],...
					'ch_ylabel','Distribution MU Rx Power');
				m_multiplot(2,3,2) = GFigure('m_X',	[min_rxPower,min_rxPower],...
					'm_Y',[0,1.2*max(v_pdfPower)]);									
				F(1) = GFigure('m_multiplot',m_multiplot);
				F.plot();
				v_frames(ind_iter) = getframe(gcf);
				%pause();
				pause(.005)
				
				
				% Gradient
				v_positionsMUMeasurement = m_positionsMU( : , randperm( num_MU , miniBatchSize ) );				
				[m_avgGradient] = UAVPlacementExperiments.stochasticGradient2D( ...
					v_txPowerBSs , m_positionsBSs , f_gain , f_gainDXBS, f_utilityPrime , v_positionsMUMeasurement );				
								
				% Gradient step
				m_positionsBSs = m_positionsBSs + m_avgGradient*diag(v_stepSizes);
				
			end
			
			
			% Create video
			writerObj = VideoWriter('./Experiments/UAVPlacementExperiments_data/aerialBaseStations',...
				'Motion JPEG AVI');
		%	'Motion JPEG 2000');
			writerObj.FrameRate = 10;
			writerObj.open();
			for ind_frame = 1:length(v_frames)
				writerObj.writeVideo(v_frames(ind_frame));				
			end
			writerObj.close();
			
		end	
	
		% Aerial BSs seek their optimum location using stochastic gradient
		% descent with minibatches.
		% - Utility: total power
		% - finitely many MUs. All same throughput
		% - relative to 2001, this experiment is restructured using a
		% simulator. 
		% Observations: the BSs nicely spread across space.
		function F = experiment_2003(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			
			regionSideLength = 7; % region is a square with bottom left edge at the origin. This is the side length
			
			% Mobile users
			num_MU = 200;		
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','totalPower',...
				'powerStart',2,'powerStop',3 );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
				
			% Base stations
			height = 1;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',15,...
				'v_position',[0,0,height]','b_fixedHeight',1);
			v_txPowerBSs = [2,3,3,3,4,4,6];
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace');
			
			% Plotter
			pp = PathPlotter('regionSideLength',regionSideLength);
			
			% Simulator
			sim = PlacementSimulator('num_steps',100,'miniBatchSize',200,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			sim.generateSamplePath();
						
			F = [];
		end	
		
		
		% Aerial BSs seek their optimum location using stochastic gradient
		% descent with minibatches.
		% - Utility: max power
		% - finitely many MUs. All same throughput
		% - similar to 2003, but max power metric rather than sum power
		% Observations: the BSs nicely spread across space.
		function F = experiment_2004(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			
			regionSideLength = 7; % region is a square with bottom left edge at the origin. This is the side length
			
			% Mobile users
			num_MU = 200;		
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',2,'powerStop',3 );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
				
			% Base stations
			height = 1;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',100,...
				'v_position',[0,0,height]','b_fixedHeight',1);
			v_txPowerBSs = 2.5*[2,3,3,3,4,4,6];
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace');
			
			% Plotter
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric','maxPower');
			
			% Simulator
			sim = PlacementSimulator('num_steps',200,'miniBatchSize',50,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			sim.generateSamplePath();
						
			F = [];
		end	
		
		% UNDER CONSTRUCTION
		% Aerial BSs seek their optimum location using stochastic gradient
		% descent with minibatches.
		% - Utility: max power
		% - finitely many MUs. All same throughput
		% - like 2004, but the numbers are more realistic --> ref. 
		% conversation with Alberto [ Maximum Coupling loss en LTE ] 
		% Observations: the BSs nicely spread across space.
		function F = experiment_2005(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			
			regionSideLength = 7e3; % region is a square with bottom left 
			%    edge at the origin. This is the side length in meters
			
			% Mobile users
			num_MU = 1;% 200;				
			powerStart = dBm_to_Watt(-91);
			powerStop = dBm_to_Watt(-80);
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'b_dBm',1,'powerStart',powerStart,'powerStop',powerStop);
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
				
			% Base stations
			height = 50;
			num_BS = 1;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',1e8,...
				'v_position',[0,0,height]','b_fixedHeight',1);
			v_txPowerBSs = dBm_to_Watt(7)*ones(num_BS,1);%250*[2,3,3,3,4,4,6];
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			f=3605e6; lambda = 3e8/f; d = 1; 
			gain1mDB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi); % in dB
			gain1m = 10^(gain1mDB/10)
			ch = Channel('ch_model','freeSpace','gainConstant',gain1m);
			
			% Plotter
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric','maxPower',...
				'b_dBm',1);
			
			% Simulator
			sim = PlacementSimulator('num_steps',200,'miniBatchSize',1,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			sim.generateSamplePath();
						
			F = [];
		end	
		
		% rescaling 2004 to get more appropriate units
		function F = experiment_2006(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			
			rng(0)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 3.971641173621411e-13; % scaling to avoid numerical problems. 
			
			% Mobile users
			num_MU = 200;
			powerStart = 2;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = 3;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',powerStart,'powerStop',powerStop );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
				
			% Base stations
			height = 0.03;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',25,...
				'v_position',[0,0,height]','b_fixedHeight',1);
			gainConstant = 10^(-94/10); % gain at 1 km
			%gainConstant = 3;
			v_txPowerBSs = 2.5*[2,3,3,3,4,4,6]/gainConstant;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',1);
			
			% Simulator
			sim = PlacementSimulator('num_steps',200,'miniBatchSize',50,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			sim.generateSamplePath();
						
			F = [];
		end	
				
		% slight modification of 2006 to get more round numbers
		function F = experiment_2007(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			
			rng(0)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 4e-13; %3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 200;
			%powerStart = 2;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',powerStart,'powerStop',powerStop );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
			m_positionsRemoteMUs = 100*regionSideLength*...
				[ 1 1 0; -1 1 0; 1 -1 0]';
			v_mobileUsers = [v_mobileUsers;...
				templateMobileUser.cloneAtPositions(m_positionsRemoteMUs)];
			
			
			% Base stations
			height = 0.03;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',10,... %25
				'v_position',[0,0,height]','b_fixedHeight',1);
			%gainConstant = 10^(-94/10); % gain at 1 km
			f=2385e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);
			
			v_txPowerBSs = dBm_to_Watt([7,9,9,9,10,10,12])/wattPerPowerUnit;
			%v_txPowerBSs = dBm_to_Watt([7,9])/wattPerPowerUnit;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			v_powerLimits =  [ dBm_to_Watt(-100)/wattPerPowerUnit,
				dBm_to_Watt(-80)/wattPerPowerUnit];
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',1,'ch_distanceUnits','km','v_powerLimitsColor',v_powerLimits,...
				'v_powerLimitsHistogram',v_powerLimits);
			
			% Simulator
			sim = PlacementSimulator('num_steps',50,'miniBatchSize',50,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			[~,F] = sim.generateSamplePath();
						
		end	
		
		% Equal to 2007 but slower (smaller step size and greater number of
		% iterations).
		function F = experiment_2008(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			displaySettings.v_windowPosition = [   744         449        1094         601];
			
			rng(0)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 4e-13; %3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 200;
			%powerStart = 2;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',powerStart,'powerStop',powerStop );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
				
			% Base stations
			height = 0.03;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',2,... %25
				'v_position',[0,0,height]','b_fixedHeight',1);
			%gainConstant = 10^(-94/10); % gain at 1 km
			f=2390e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);
			
			v_txPowerBSs = dBm_to_Watt([7,9,9,9,10,10,12])/wattPerPowerUnit;
			%v_txPowerBSs = dBm_to_Watt([7,9])/wattPerPowerUnit;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			v_powerLimits =  [ dBm_to_Watt(-100)/wattPerPowerUnit,
				dBm_to_Watt(-80)/wattPerPowerUnit];
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',1,'ch_distanceUnits','km','v_powerLimitsColor',v_powerLimits,...
				'v_powerLimitsHistogram',v_powerLimits);
			
			% Simulator
			sim = PlacementSimulator('num_steps',250,'miniBatchSize',50,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			[~,F] = sim.generateSamplePath();
						
		end	
			
		% Experiment with kmeans.
		function F = experiment_2009(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			
			rng(0)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 4e-13; %3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 200;
			%powerStart = 2;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','kmeans',...
				'powerStart',powerStart,'powerStop',powerStop );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
			
			m_positionsRemoteMUs = 100*regionSideLength*...
				[ 1 1 0; -1 1 0; 1 -1 0]';
			v_mobileUsers = [v_mobileUsers;...
				templateMobileUser.cloneAtPositions(m_positionsRemoteMUs)];
			
			% Base stations
			height = 0.03;
			templateAerialBaseStation = KmeansAerialBaseStation('stepSize',2,... %25
				'v_position',[0,0,height]','b_fixedHeight',1);
			%gainConstant = 10^(-94/10); % gain at 1 km
			f=2390e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);
			
			%v_txPowerBSs = dBm_to_Watt([7,9,9,9,10,10,12])/wattPerPowerUnit;
			v_txPowerBSs = dBm_to_Watt([3,9,3,9,10,10,13])/wattPerPowerUnit;
			%v_txPowerBSs = dBm_to_Watt([7,9])/wattPerPowerUnit;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			v_powerLimits =  [ dBm_to_Watt(-100)/wattPerPowerUnit,
				dBm_to_Watt(-80)/wattPerPowerUnit];
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',1,'ch_distanceUnits','km','v_powerLimitsColor',v_powerLimits,...
				'v_powerLimitsHistogram',v_powerLimits);
			
			% Simulator
			sim = PlacementSimulator('num_steps',250,'miniBatchSize',num_MU,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			[~,F] = sim.generateSamplePath();
						
		end	
		
		% Comparison of histograms: Stochastic vs. kmeans
		% Experiment for Globecom
		function F = experiment_2010(obj,niter)
				
			global displaySettings 
			displaySettings.v_windowPosition = [   744         449        1094         601];
			
			rng(0)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 4e-13; %3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 200;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',powerStart,'powerStop',powerStop );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
			m_positionsRemoteMUs = 5*regionSideLength*...
				[ 1 1 0; -1 1 0]';
			%	[ 1 1 0; -1 1 0; 1 -1 0]';
			v_mobileUsers = [v_mobileUsers;...
				templateMobileUser.cloneAtPositions(m_positionsRemoteMUs)];
			
			% Base stations
			height = 0.03;			
			f=2385e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);			
			v_txPowerBSs = dBm_to_Watt([7,9,9,9,12])/wattPerPowerUnit;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)			
			% 1. Stochastic BSs
			rng(1); % to ensure same intial positions for both sets of AirBSs
			templateStochAerialBaseStation = StochasticAerialBaseStation('stepSize',5,... 
				'v_position',[0,0,height]','b_fixedHeight',1);
			v_stochAerialBaseStations = templateStochAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);			
			% 2. Kmeans BSs
			rng(1); % to ensure same intial positions for both sets of AirBSs
			templateKmeansAerialBaseStation = KmeansAerialBaseStation(... 
				'v_position',[0,0,height]','b_fixedHeight',1);
			v_kmeansAerialBaseStations = templateKmeansAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
			
			% Channel
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			v_powerLimits =  [ dBm_to_Watt(-100)/wattPerPowerUnit,
				dBm_to_Watt(-80)/wattPerPowerUnit];
			ppTemplate = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',1,'ch_distanceUnits','km','v_powerLimitsColor',v_powerLimits,...
				'v_powerLimitsHistogram',v_powerLimits);	
			%ppTemplate.ch_videoFileName = '';
			ppStoch = clone(ppTemplate);
			ppKmeans = clone(ppTemplate);
			
			% Simulator
			sim = PlacementSimulator('num_steps',100,'miniBatchSize',50,...  % number of measurements per update;
				'channel',ch);
			
			% Simulation		
			% 1. Stochastic AirBSs
			sim.v_mobileUsers = v_mobileUsers;
			sim.v_aerialBaseStations = v_stochAerialBaseStations;
			sim.pathPlotter = ppStoch;
			rng(2) % ensure same simulation setup
			[~,Fstoch] = sim.generateSamplePath();
				
			% 2. Kmeans AirBSs
			for ind_MU = 1:length(v_mobileUsers)
				v_mobileUsers(ind_MU).ch_utility = 'kmeans';
			end
			sim.v_mobileUsers = v_mobileUsers;
			sim.v_aerialBaseStations = v_kmeansAerialBaseStations;
			sim.pathPlotter = ppKmeans;
			rng(2)
			[~,Fkmeans] = sim.generateSamplePath();
			
			Fcompare = ppStoch.compareHistograms( ppKmeans );
			
			F = [Fstoch(1),Fkmeans(1),Fcompare];
			
		end	
		
		% Video for youtube with different number of BSs
		function F = experiment_2011(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			displaySettings.v_windowPosition = [   44         449        1094         601];
			
			rng(0)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 4e-13; %3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 200;
			%powerStart = 2;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',powerStart,'powerStop',powerStop );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1);
				
			% Base stations
			height = 0.03;
			templateAerialBaseStation = StochasticAerialBaseStation(... %25
				'v_position',[0,0,height]','b_fixedHeight',1);
			%gainConstant = 10^(-94/10); % gain at 1 km
			f=2390e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);
			
			c_v_txPowerBSs = {dBm_to_Watt([7,9,9,9,12])/wattPerPowerUnit,...
				dBm_to_Watt([7,10,10,12,12])/wattPerPowerUnit,...
				dBm_to_Watt(fliplr([7,7,7,7,10,10,12]))/wattPerPowerUnit,...
				dBm_to_Watt([16,16])/wattPerPowerUnit,...
				dBm_to_Watt(fliplr([7,7,7,7,7,8,8,8,8,8]-1))/wattPerPowerUnit};
			
			% Channel
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			v_powerLimits =  [ dBm_to_Watt(-100)/wattPerPowerUnit,
				dBm_to_Watt(-80)/wattPerPowerUnit];
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',1,'ch_distanceUnits','km','v_powerLimitsColor',v_powerLimits,...
				'v_powerLimitsHistogram',v_powerLimits,'frameRate',20);
			%pp.ch_videoFileName = '';
			
			% Simulator
			sim = PlacementSimulator('miniBatchSize',50,...
				'v_mobileUsers',v_mobileUsers,...
				'channel',ch,'pathPlotter',pp); % number of measurements per update;
				
			% Simulation
			v_frames = [];
			v_numSteps = [200,200,250,100,500];
			%v_numSteps = [2,2,2,1,400];
			v_stepSize = [2,  3,  3,  3,  6  ];
			for ind_setting = 1:length(c_v_txPowerBSs)
				v_txPowerBSs = c_v_txPowerBSs{ind_setting};
				v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit);
				templateAerialBaseStation.stepSize = v_stepSize(ind_setting);
				sim.v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1);
				sim.num_steps = v_numSteps(ind_setting);
				
				[~,v_F] = sim.generateSamplePath();
				v_F(1).plot()
				v_frames = [v_frames, pp.v_frames];
				%pause()
			end
			
			if pp.ch_videoFileName
				pp.v_frames = v_frames;
				pp.writeVideo();
				
				load handel
				sound(y(1:17000),Fs)
			end
			F = [];
		end	
		
		% Experiment with leaky ReLU
		function F = experiment_2012(obj,niter)
				
			global displaySettings 
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			%displaySettings.v_windowPosition = [   744         449        1094         601];
			displaySettings.v_windowPosition = [   744-500        449        1094         601];
			
			rng(1)
			
			regionSideLength = 7; % km % region is a square with bottom left edge at the origin. This is the side length
			wattPerPowerUnit = 1; %3.971641173621411e-13; % scaling to 
			%wattPerPowerUnit = 3.971641173621411e-13; % scaling to 
			%    avoid numerical problems. One natural power unit henceforth equals $wattPerPowerUnit$ Watt.
			
			% Mobile users
			num_MU = 2;
			b_2D = 0;
			%powerStart = 2;
			powerStart = dBm_to_Watt(-91)/wattPerPowerUnit;
			powerStart_dBm = Watt_to_dBm( powerStart*wattPerPowerUnit )
			powerStop = dBm_to_Watt(-89)/wattPerPowerUnit;
			powerStop_dBm = Watt_to_dBm( powerStop*wattPerPowerUnit )
			% templateMobileUser = MobileUser('ch_utility','leakyRelu','ch_aggregate','smoothMaxPower',...
			% 	'powerStart',powerStart , 'b_dBm',1 );
			templateMobileUser = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'powerStart',powerStart ,'powerStop',powerStop, 'b_dBm',1, 'positiveScaling',0.0 );
			v_mobileUsers = templateMobileUser.cloneAtRandomPositions(num_MU,regionSideLength,1,b_2D);
				
			% Base stations
			height = 0.03;
			templateAerialBaseStation = StochasticAerialBaseStation('stepSize',1e4,... %25
				'v_position',[0,0,height]','b_fixedHeight',1,'b_debug',1);
			%v_txPowerBSs = dBm_to_Watt([7,9,9,9,10,10,12]-3)/wattPerPowerUnit;
			v_txPowerBSs = dBm_to_Watt([4,6])/wattPerPowerUnit;
			v_txPowerBSs_dBm = Watt_to_dBm(v_txPowerBSs*wattPerPowerUnit)
			v_aerialBaseStations = templateAerialBaseStation.CloneAtRandomPositions(v_txPowerBSs,regionSideLength/3,1,b_2D);
			
			% Channel
			f=2390e6; lambda = 3e8/f/1000; d = 1;  % lambda and d in km
			gainConstant_dB = 6 ... % AirBS only radiates downwards
				+ 0 ... % isotropic UE antenna
				- 20*log10( d/lambda) - 20*log10(4*pi) % gain @1km in dB
			gainConstant = 10^(gainConstant_dB/10);
			ch = Channel('ch_model','freeSpace','gainConstant',gainConstant);
			
			% Plotter
			%v_powerLimits =  [ dBm_to_Watt(-100)/wattPerPowerUnit,
			% dBm_to_Watt(-80)/wattPerPowerUnit];
			v_powerLimits = [powerStart,powerStop];
			pp = PathPlotter('regionSideLength',regionSideLength,'ch_metric',...
				'maxPower','wattPerPowerUnit', wattPerPowerUnit,...
				'b_dBm',0,'ch_distanceUnits','km','v_powerLimitsColor',v_powerLimits,...
				'v_powerLimitsHistogram',v_powerLimits);
			
			% Simulator
			sim = PlacementSimulator('num_steps',250,'miniBatchSize',NaN,...
				'v_mobileUsers',v_mobileUsers,'v_aerialBaseStations',v_aerialBaseStations,...
				'channel',ch,'pathPlotter',pp,'b_debug',1); % number of measurements per update;
							
			% Simulation
			[~,F] = sim.generateSamplePath();
						
		end	
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 99. Additional auxiliary experiments
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% Plot of sigmoid and shifted-scaled sigmoid.
		function F = experiment_9901(obj,niter)
			
			global displaySettings
			displaySettings.v_windowPosition = [100   402   444   123];
			
			v_x = linspace(-5,5,1000);			
			
			% Functions
			v_unitStep = (v_x > 0);
			v_sigmoid = sigmoid(v_x);
			
			mu = MobileUser('powerStart',0,'powerStop',1);
			v_modifiedSigmoid = mu.shiftedSigmoid( v_x );
		
			F = GFigure('m_X',v_x,'m_Y',[v_unitStep;v_sigmoid;v_modifiedSigmoid],...
				'ch_xlabel','Abscissa, x','ch_ylabel','Function value',...
				'c_legend',{'UnitStep(x)','Sigmoid(x)','ModifiedSigmoid(x)'},...
				'c_styles',{'-','--','-.'},'ch_legendPosition','NorthWest');
		end
		
		% Plot of MU utility and its derivative in dBm
		function F = experiment_9902(obj,niter)
			
			v_power_dBm = linspace(-200,-50,1000); 	 
			
			powerStart = dBm_to_Watt(-91);
			powerStop = dBm_to_Watt(-80);
			mu = MobileUser('ch_utility','sigmoid','ch_aggregate','smoothMaxPower',...
				'b_dBm',1,'powerStart',powerStart,'powerStop',powerStop);
			for ind = length(v_power_dBm):-1:1
				v_utility(ind) = mu.utility(dBm_to_Watt(v_power_dBm(ind)));
				v_utilityGradient(ind) = 0*mu.utilityGradient(dBm_to_Watt(v_power_dBm(ind)));
			end
		
			F = GFigure('m_X',dBm_to_Watt(v_power_dBm),'m_Y',[v_utility;v_utilityGradient],...
				'ch_xlabel','Abscissa, x [dBm]','ch_ylabel','Function value',...		
				'c_legend',{'UnitStep(x)','Sigmoid(x)','ModifiedSigmoid(x)'},...
				'c_styles',{'-','--','-.'},'ch_legendPosition','NorthWest');
		end
		
		
	end
		
	methods(Static) % OLD FUNCTIONS
		
		function [v_totalPower1D , m_partialPower] = totalPower1D( v_txPowerBSs , v_positionsBS , f_gain , v_positionsMU )
			% v_positionsBS and v_positionsMU
			assert(size(v_positionsMU,1) == 1 ); % row vector
			num_MU = length(v_positionsMU);
			num_BS = length(v_positionsBS);
			m_partialPower = diag(v_txPowerBSs)*f_gain(repmat(v_positionsBS,1,num_MU),...
				repmat(v_positionsMU,num_BS,1));
			v_totalPower1D = sum(m_partialPower,1);
		end
		
		function [v_totalPower , m_partialPower] = totalPower2D( v_txPowerBSs , m_positionsBS , f_gain , m_positionsMU )
			% INPUT:
			% m_positionsBS : 2 x num_BS
			% m_positionsMU : 2 x num_MU
			% f_gain takes 2 arguments with 2 rows each
			%
			% OUTPUT:
			% m_totalPower : 1 x num_MU matrix
			% m_partialPower : num_BS x num_MU matrix
			%
			assert(size(m_positionsBS,1) == 2 ); 
			assert(size(m_positionsMU,1) == 2 ); 
			num_BS = size(m_positionsBS,2);
			for ind_BS = num_BS:-1:1
				m_partialPower(ind_BS,:) = v_txPowerBSs(ind_BS)*f_gain(m_positionsBS(:,ind_BS),...
				m_positionsMU);
			end
			v_totalPower = sum(m_partialPower,1);
		end
		
		function [v_avgGradient, m_instantaneousGradient] = stochasticGradient1D( v_txPowerBSs , v_positionsBS , f_gain , f_gainDXBS, f_utilityPrime , v_positionsMUMeasurement )
			% v_avgGradient is num_BS x 1
			% m_instantaneousGradient is num_BS x num_MU, where num_MU = length(v_positionsMUMeasurement);
			
			assert(size(v_positionsMUMeasurement,1) == 1 ); % row vector
			num_MU = length(v_positionsMUMeasurement);
			num_BS = length(v_positionsBS);
			
			[v_totalPower1D] = UAVPlacementExperiments.totalPower1D( v_txPowerBSs , v_positionsBS , f_gain , v_positionsMUMeasurement );
			m_instantaneousGradient = repmat(f_utilityPrime( v_totalPower1D ),num_BS,1)...
					.* (diag(v_txPowerBSs)* f_gainDXBS( repmat(v_positionsBS,1,num_MU) , ...
					repmat(v_positionsMUMeasurement,num_BS,1)));%      one row per BS, one col per MU
			v_avgGradient = mean(m_instantaneousGradient,2);
			
		end
			
		function [m_avgGradient, t_instantaneousGradient] = stochasticGradient2D( v_txPowerBSs , m_positionsBS , f_gain , f_gainDXBS, f_utilityPrime , m_positionsMUMeasurement )
			% INPUT:
			% m_positionsBS : 2 x num_BS
			% m_positionsMUMeasurement : 2 x num_MU
			% f_gain takes 2 arguments with 2 rows each
			%
			% OUTPUT:
			% v_avgGradient is 2 x num_BS
			% t_instantaneousGradient is 2 x num_BS x num_MU
			%
			assert(size(m_positionsMUMeasurement,1) == 2 );
			assert(size(m_positionsBS,1) == 2 );
			num_BS = size(m_positionsBS,2);
			
			[m_totalPower] = UAVPlacementExperiments.totalPower2D( v_txPowerBSs , m_positionsBS , f_gain , m_positionsMUMeasurement );
			for ind_BS = num_BS:-1:1
				v_positionBSNow = m_positionsBS(:,ind_BS);
				m_gradientsThisBS =  f_gainDXBS( v_positionBSNow , m_positionsMUMeasurement)...
					*diag(f_utilityPrime( m_totalPower )); % one row per dimension, one column per MU
				
				t_instantaneousGradient(:,ind_BS,:) = permute(m_gradientsThisBS,[1 3 2]);
			end
			m_avgGradient = mean(t_instantaneousGradient,3);
			
		end
		
		
	end 
	
	
end
