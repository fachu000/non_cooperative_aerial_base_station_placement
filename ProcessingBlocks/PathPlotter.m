classdef PathPlotter < ProcessingBlock
	%PathPlotter Summary of this class goes here
	%   Detailed explanation goes here
	

	properties
			
		regionSideLength
		v_plotMobileUsers = MobileUser(); % Test mobile users at each grid point
		v_fractionUsersOverMinMetric = [];
		
		m_xCoordinatesMUPlot,m_yCoordinatesMUPlot
		
		% Figure selection		
		ch_metric = 'totalPower'; % Possible values:
		%    - 'totalPower'
		%    - 'maxPower'
		b_dBm = 0; % set to 1 for plots in dBm
		
		num_pointsGridPerDim = 30;
		ch_videoFileName = './Experiments/UAVPlacementExperiments_data/aerialBaseStations';
		%    Set to ch_videoFileName = [] to skip the video creation
		frameRate = 10;		
		v_frames = [];
		
		% Temporal solution: power units are scaled for representation
		%     purposes
		wattPerPowerUnit = 1;
		ch_distanceUnits = 'm'; % can be e.g. 'km'
			
		m_path = []; % stored path of base stations
		
		v_powerLimitsColor = []; % if set to a vector of 2 parameters, then this 
		%    will be the limits of the power for the figures. 
		v_powerLimitsHistogram = [];
		v_powerLimitsSurf = [];
		
		% initial histogram, of the form [v_binCenters;v_num_MUEachBin] (2
		% rows)
		m_initialHistogram = [];
		m_finalHistogram = []; % idem 
		
	end
	
	methods
		
		% Constructor
		function obj = PathPlotter(varargin)
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end
		end
				
		function obj = initialize(obj,templateTestMobileUser)
			% templateTestMobileUser :object of class MobileUser. It will
			%     be replicated at all grid points.
			
			assert(numel(templateTestMobileUser)==1);
			obj.m_path = [];
			obj.m_initialHistogram = [];
			obj.m_finalHistogram = [];
			obj.v_frames = [];
			obj.v_fractionUsersOverMinMetric = [];
			obj.v_powerLimitsSurf = [];
			
			% Preparations for plots			
			[obj.m_xCoordinatesMUPlot,obj.m_yCoordinatesMUPlot] = meshgrid(...
				linspace(0,obj.regionSideLength,obj.num_pointsGridPerDim),...
				linspace(0,obj.regionSideLength,obj.num_pointsGridPerDim)); 
			%m_yCoordinatesMUPlot = flipud(m_yCoordinatesMUPlot); % Y axis increases upwards
			m_positionsPlotMUs = [obj.m_xCoordinatesMUPlot(:),obj.m_yCoordinatesMUPlot(:),zeros(size(obj.m_yCoordinatesMUPlot(:)))];
			
			for ind_plotMU = 1:size(m_positionsPlotMUs,1)
				obj.v_plotMobileUsers(ind_plotMU,1) = clone(templateTestMobileUser);
				obj.v_plotMobileUsers(ind_plotMU,1).v_position = m_positionsPlotMUs(ind_plotMU,:)';				
			end
			
			
		end
			
		function obj = receiveNewPlacement(obj,v_aerialBaseStations,channel,v_mobileUsers)
			
			[m_positionsBSs,m_positionsMUs,ch_metricName,v_metricPlot,...
				v_metric,v_utility,v_num_MUEachBin,v_binCenters,...
				powerStart,powerStop,ch_xCoordinate,ch_yCoordinate,v_caxisPower] ...
				= collectMetricsToPlot(obj,v_aerialBaseStations,channel,v_mobileUsers);			
			
			if isempty(obj.ch_videoFileName)
				return
			end
			
			%% 1. Plot of AirBS and MU positions
			m_multiplot(1,1,1) = GFigure('m_X',m_positionsMUs(1,:),'m_Y',m_positionsMUs(2,:),...
				'm_Z',m_positionsMUs(3,:),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_zlabel',['Height [',obj.ch_distanceUnits,']'],...
				'ch_plotType3D','stem3','c_styles',{'.k'},...
				'ch_title','Locations of mobile users and aerial base stations');
			%'c_legend',{'Mobile Users','Aerial BSs'});
			m_multiplot(1,1,2) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
				'm_Z',m_positionsBSs(3,:),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_plotType3D','stem3','c_styles',{'sb'});
			
			%% 2. 3D view (surf) of metric	
			if isempty(obj.v_powerLimitsSurf)
				obj.v_powerLimitsSurf = [0.9*min(min(v_metricPlot)),1.2*max(max(v_metricPlot))];
			end
			m_multiplot(1,2,1) = GFigure('m_X',obj.m_xCoordinatesMUPlot,...
				'm_Y',obj.m_yCoordinatesMUPlot,...
				'm_Z',reshape(v_metricPlot,obj.num_pointsGridPerDim,obj.num_pointsGridPerDim),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_zlabel',ch_metricName,...
				'ch_plotType3D','surf','v_zlim',obj.v_powerLimitsSurf,...
				'v_caxis',obj.v_powerLimitsSurf);
			
			%% 3. Fraction of MUs with metric over powerStart
			ch_title = sprintf('No. served users (# users over Pmin): %d / %d',round(obj.v_fractionUsersOverMinMetric(end)*length(v_mobileUsers)),length(v_mobileUsers));
			m_multiplot(1,3,1) = GFigure('m_X',1:length(obj.v_fractionUsersOverMinMetric),...
				'm_Y',[obj.v_fractionUsersOverMinMetric],...
				'ch_xlabel','Iteration index',...
				'c_legend',{'Frac. area over Pmin'},'ch_legendPosition','northwest',...
				'v_ylim',[0,1],'ch_title',...
				ch_title);
			
			%% 4. Color map of metric			
			m_multiplot(2,1,1) = GFigure('m_X',obj.m_xCoordinatesMUPlot,...
				'm_Y',obj.m_yCoordinatesMUPlot,...
				'm_Z',reshape(v_metricPlot,obj.num_pointsGridPerDim,obj.num_pointsGridPerDim),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...				
				'ch_plotType3D','imagesc','v_caxis',v_caxisPower,...
				'ch_title',ch_metricName);
			m_multiplot(2,1,2) = GFigure('m_X',m_positionsMUs(1,:),'m_Y',m_positionsMUs(2,:),...
				'm_Z',m_positionsMUs(3,:),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_plotType3D','stem3','c_styles',{'.w'});
			m_multiplot(2,1,3) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
				'm_Z',m_positionsBSs(3,:),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_plotType3D','stem3','c_styles',{'sc'},'m_colorset',[1,1,1]);
			
			%% 5. 3D stem of the utility at each MU		
			ch_title = sprintf('Avg. utility = %g',mean(v_utility));
			m_multiplot(2,2) = GFigure('m_X',m_positionsMUs(1,:),'m_Y',m_positionsMUs(2,:),...
				'm_Z',v_utility,...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_zlabel','Utility',...
				'ch_title',ch_title,...
				'ch_zlabel','Utility','ch_plotType3D','stem3'); %,'c_styles','filled');
			
			%% 6. Histogram of metric
			v_pdfMetric = v_num_MUEachBin; % /length(v_metric)/(v_binCenters(2)-v_binCenters(1));
			m_multiplot(2,3,1) = GFigure('m_X',v_binCenters,...
				'm_Y',v_pdfMetric,...
				'c_legend',{'Distribution','Pmin'},'ch_legendPosition','northwest',...
				'ch_xlabel',ch_metricName,'v_ylim',[0,1.4*max(obj.m_initialHistogram(2,:))],...
				'ch_ylabel','Histogram [no. mobile users]','ch_plotType2D','stem');
			m_multiplot(2,3,2) = GFigure('m_X',	[powerStart,powerStart],...
				'm_Y',[0,1.4*max(obj.m_initialHistogram(2,:))],'m_colorset',[0,0.5,0]);
						
			%% Plot and store
			F(1) = GFigure('m_multiplot',m_multiplot);
			F.plot();
			obj.v_frames = [obj.v_frames,getframe(gcf)];
			%pause();
			pause(.005)
			
			
		end
		
		function F = finalize(obj,v_aerialBaseStations,channel,v_mobileUsers)
			
			% Create video
			if ~isempty(obj.ch_videoFileName)
				obj.writeVideo();
			end
			
			% Power Map figure
			[m_positionsBSs,m_positionsMUs,ch_metricName,v_metricPlot,...
				v_metric,v_utility,v_num_MUEachBin,v_binCenters,...
				powerStart,powerStop,ch_xCoordinate,ch_yCoordinate,v_caxisPower] ...
				= collectMetricsToPlot(obj,v_aerialBaseStations,channel,v_mobileUsers);
			obj.m_finalHistogram = [v_binCenters;v_num_MUEachBin];
			
			m_multiplot(1,1,1) = GFigure('m_X',obj.m_xCoordinatesMUPlot,...
				'm_Y',obj.m_yCoordinatesMUPlot,...
				'm_Z',reshape(v_metricPlot,obj.num_pointsGridPerDim,obj.num_pointsGridPerDim),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...				
				'ch_plotType3D','imagesc','v_caxis',v_caxisPower,...
				'ch_title',['Clipped ',ch_metricName],'b_colorbar',1);
			m_multiplot(1,1,2) = GFigure('m_X',permute(obj.m_path(1,:,:),[2 3 1]),'m_Y',permute(obj.m_path(2,:,:),[2 3 1]),...
				'm_Z',permute(obj.m_path(3,:,:),[2 3 1]),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_plotType3D','plot3','c_styles',{'-'},...
				'm_colorset',[0,0.5,0]);
			m_multiplot(1,1,3) = GFigure('m_X',m_positionsMUs(1,:),'m_Y',m_positionsMUs(2,:),...
				'm_Z',m_positionsMUs(3,:),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_plotType3D','stem3','c_styles',{'.w'},...
				'm_colorset',[0,0,0.5]);
			m_multiplot(1,1,4) = GFigure('m_X',m_positionsBSs(1,:),'m_Y',m_positionsBSs(2,:),...
				'm_Z',m_positionsBSs(3,:),...
				'ch_xlabel',ch_xCoordinate,'ch_ylabel',ch_yCoordinate,...
				'ch_plotType3D','stem3','c_styles',{'sc'});
			F(1) = GFigure('m_multiplot',m_multiplot);
			
			%% Histogram of metric
			m_multiplot = GFigure();
			v_binCenters_ini = [obj.m_initialHistogram(1,:)];
			v_num_MUEachBin_ini = [obj.m_initialHistogram(2,:)];
			
			%v_pdfMetric = v_num_MUEachBin/length(v_metric)/(v_binCenters(2)-v_binCenters(1));
			m_multiplot(1,1,1) = GFigure('m_X',v_binCenters_ini,...
				'm_Y',v_num_MUEachBin_ini,...
				'c_legend',{'Distribution','Pmin'},...
				'ch_xlabel',ch_metricName,'v_ylim',[0,1.2*max(v_num_MUEachBin_ini)],...
				'ch_ylabel','Histogram (no. MUs)','ch_plotType2D','stem',...
				'ch_title','Initial histogram');
			m_multiplot(1,1,2) = GFigure('m_X',	[powerStart,powerStart],...
				'm_Y',[0,1.2*max(v_num_MUEachBin_ini)],'m_colorset',[0,0.5,0]);
			
			m_multiplot(2,1,1) = GFigure('m_X',v_binCenters,...
				'm_Y',v_num_MUEachBin,...
				'c_legend',{'Distribution','Pmin'},...
				'ch_xlabel',ch_metricName,'v_ylim',[0,1.2*max(v_num_MUEachBin)],...
				'ch_ylabel','Histogram (no. MUs)','ch_plotType2D','stem',...
				'ch_title','Final histogram');
			m_multiplot(2,1,2) = GFigure('m_X',	[powerStart,powerStart],...
				'm_Y',[0,1.2*max(v_num_MUEachBin)],'m_colorset',[0,0.5,0]);
				
			F(2) = GFigure('m_multiplot',m_multiplot);
			
		end
			
		function writeVideo(obj)
			
			writerObj = VideoWriter(obj.ch_videoFileName,...
				'Motion JPEG AVI');
			%	'Motion JPEG 2000');
			writerObj.FrameRate = obj.frameRate;
			writerObj.open();
			for ind_frame = 1:length(obj.v_frames)
				writerObj.writeVideo(obj.v_frames(ind_frame));
			end
			writerObj.close();
			fprintf('Video stored in file %s\n',obj.ch_videoFileName);
		end
		
		function powerStart = getPowerStart(obj)
			
			powerStart = obj.v_plotMobileUsers(1,1).powerStart*obj.wattPerPowerUnit;
			
			if obj.b_dBm
				powerStart = Watt_to_dBm(powerStart);
			end
			
		end
		
		function F = compareHistograms( obj , benchmarkPathPlotter )
			
			powerStart = getPowerStart(obj);
			ch_metricName = getMetricName(obj);
			
			%% Histogram of metric
			m_multiplot = GFigure();
			v_binCenters_ini = obj.m_initialHistogram(1,:);
			v_num_MUEachBin_ini = obj.m_initialHistogram(2,:);
			
			v_binCenters1 = obj.m_finalHistogram(1,:);
			v_num_MUEachBin1 = obj.m_finalHistogram(2,:);
			fprintf('Number of users not covered by the proposed method: %d\n',sum( v_num_MUEachBin1(v_binCenters1<powerStart) ));
			v_binCenters2 = benchmarkPathPlotter.m_finalHistogram(1,:);
			v_num_MUEachBin2 = benchmarkPathPlotter.m_finalHistogram(2,:);
			fprintf('Number of users not covered by the benchmark method: %d\n',sum( v_num_MUEachBin2(v_binCenters2<powerStart) ));
			
			v_ylim = [0,1.2*max([v_num_MUEachBin1, v_num_MUEachBin2])];
			
			%v_pdfMetric = v_num_MUEachBin/length(v_metric)/(v_binCenters(2)-v_binCenters(1));
			m_multiplot(1,1,1) = GFigure('m_X',v_binCenters_ini,...
				'm_Y',v_num_MUEachBin_ini,...
				'c_legend',{'Distribution','Pmin'},...
				'ch_xlabel',ch_metricName,'v_ylim',[0,1.2*max(v_num_MUEachBin_ini)],...
				'ch_ylabel','Histogram (no. MUs)','ch_plotType2D','stem',...
				'ch_title','Initial histogram');
			m_multiplot(1,1,2) = GFigure('m_X',	[powerStart,powerStart],...
				'm_Y',[0,1.2*max(v_num_MUEachBin_ini)],'m_colorset',[0,0.5,0]);
			
			m_multiplot(2,1,1) = GFigure('m_X',v_binCenters1,...
				'm_Y',v_num_MUEachBin1,...
				'ch_xlabel',ch_metricName,...
				'ch_ylabel','Histogram (no. MUs)','ch_plotType2D','stem',...
				'ch_title','Final histogram',...
				'c_legend',{'Proposed','Benchmark','Pmin'});
			m_multiplot(2,1,2) = GFigure('m_X',v_binCenters2,...
				'm_Y',v_num_MUEachBin2,...
				'ch_xlabel',ch_metricName,'v_ylim',v_ylim,...
				'ch_ylabel','Histogram (no. MUs)','ch_plotType2D','stem',...
				'ch_title','Final histogram','m_colorset',[0.5,0.2,0],'c_styles',{'x'});
			m_multiplot(2,1,3) = GFigure('m_X',	[powerStart,powerStart],...
				'm_Y',[0,1.2*max(v_num_MUEachBin1)],'m_colorset',[0,0.5,0]);
									
			F = GFigure('m_multiplot',m_multiplot);
			
		end
		
		
		function ch_metricName = getMetricName(obj)
			
			switch obj.ch_metric
				case 'totalPower'
					ch_metricName = 'Total Power';
				case 'maxPower'
					ch_metricName = 'Max Power';
				otherwise
					error('Not implemented')
			end
			
		end
		
	end
	
	methods(Access=private)
		
		function [m_positionsBSs,m_positionsMUs,ch_metricName,v_metricPlot,...
				v_metric,v_utility,v_num_MUEachBin,v_binCenters,...
				powerStart,powerStop,ch_xCoordinate,ch_yCoordinate,v_caxisPower] ...
				= collectMetricsToPlot(obj,v_aerialBaseStations,channel,v_mobileUsers)
			
			m_positionsBSs = v_aerialBaseStations.positionAerialBaseStations();			
			m_positionsMUs = v_mobileUsers.positionMobileUsers;
			
			% Collect metrics to plot
			[~,v_totalPowerPlot,v_maxPowerPlot] = obj.utilityAtMUs(v_aerialBaseStations,channel,obj.v_plotMobileUsers);
			[v_utility,v_totalPower,v_maxPower] = obj.utilityAtMUs(v_aerialBaseStations,channel,v_mobileUsers);
			ch_metricName = obj.getMetricName();
			switch obj.ch_metric
				case 'totalPower'					
					v_metricPlot = v_totalPowerPlot*obj.wattPerPowerUnit;
					v_metric = v_totalPower*obj.wattPerPowerUnit;
				case 'maxPower'					
					v_metricPlot = v_maxPowerPlot*obj.wattPerPowerUnit;
					v_metric = v_maxPower*obj.wattPerPowerUnit;
				otherwise
					error('Not implemented')
			end
			
			powerStart = obj.getPowerStart();
			powerStop = obj.v_plotMobileUsers(1,1).powerStop*obj.wattPerPowerUnit;
			if ~isempty(obj.v_powerLimitsColor)
				v_caxisPower = obj.v_powerLimitsColor*obj.wattPerPowerUnit;
			else
				v_caxisPower = [powerStart,powerStop]*obj.wattPerPowerUnit;					
			end
			v_powerLimitsHistogram_h = obj.v_powerLimitsHistogram*obj.wattPerPowerUnit;
			
			if obj.b_dBm
				ch_metricName = [ch_metricName ' [dBm]'];
				v_metricPlot = Watt_to_dBm( v_metricPlot);
				v_metric = Watt_to_dBm( v_metric);
				powerStop = Watt_to_dBm(powerStop);
				v_caxisPower = Watt_to_dBm(v_caxisPower);
				v_powerLimitsHistogram_h = Watt_to_dBm(v_powerLimitsHistogram_h);
				%obj.v_powerLimits = Watt_to_dBm(obj.v_powerLimits);
				%v_binCenters = linspace(min(v_metric),max(v_metric),15);
				
				if ~isempty(v_powerLimitsHistogram_h)
					v_binCenters = linspace(v_powerLimitsHistogram_h(1),v_powerLimitsHistogram_h(2),15);
					[v_num_MUEachBin] = hist( v_metric , v_binCenters);
				else					
					[v_num_MUEachBin,v_binCenters] = hist( v_metric );
				end
				%
				
			else
				ch_metricName = [ch_metricName ' [W]'];
				v_binCenters = linspace(0,max(v_metric),15);
				[v_num_MUEachBin] = hist( v_metric , v_binCenters);
			end
			
			obj.v_fractionUsersOverMinMetric = [obj.v_fractionUsersOverMinMetric,...
				mean(v_metricPlot>=powerStart)];
			ch_xCoordinate = ['x coordinate [',obj.ch_distanceUnits,']'];
			ch_yCoordinate = ['y coordinate [',obj.ch_distanceUnits,']'];
			
			if isempty(obj.m_path)
				% we are in the first iteration; let us store the histogram
				obj.m_initialHistogram = [v_binCenters;v_num_MUEachBin];
			end
			
			
			
			obj.m_path = cat(3,obj.m_path,m_positionsBSs);
		end
		
	end
	
	
	
	
	methods(Static)
		
		function [v_utility,v_totalPower,v_maxPower] = utilityAtMUs(v_aerialBaseStations,channel,v_mobileUsers)
			% v_utility, v_totalPower, and v_maxPower :  num_MU x 1 vector
			% 
			
			m_positionsBSs = v_aerialBaseStations.positionAerialBaseStations();	
			v_txPower = v_aerialBaseStations.powerAerialBaseStations();
			for ind_MU = length(v_mobileUsers):-1:1
				thisMU = v_mobileUsers(ind_MU);
				
				% Utility gradients w.r.t. power from MUs
				v_rxPower = v_txPower.*channel.gain(m_positionsBSs,thisMU.v_position);
				v_utility(ind_MU,1) = thisMU.utility(v_rxPower);
				v_totalPower(ind_MU,1) = sum(v_rxPower);
				v_maxPower(ind_MU,1) = max(v_rxPower);
				%    One row per base station, one column per MU
			end
			
		end
	
		
		
	end

end
