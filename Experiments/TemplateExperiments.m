%
%  Template for experiment files. Provide here general information on your
%  simulation study, e.g.: This file contains experiments for
%  spatio-temporal prediction of time series with drilling data. The
%  experiments 1002, 3004, 4003 and 4004 correspond to the paper ...
%
%

classdef TemplateExperiments < ExperimentFunctionSet
	
	properties
		% You may define here constants that you will use in the
		% experiments of this file
		
	end
	
	methods
				
		% Example of experiment that does not display any figure. 
		function F = experiment_1001(obj,niter)
			disp('This is a dummy experiment!!')
			A = randn(4)
			
			F = [];
		end	
		
		% Example of experiment that displays a figure. The figure is
		% first constructed and displayed through MATLAB standard commands.
		% The figure is then stored inside an GFigure object through the
		% function GFigure.captureCurrentFigure(); 
		function F = experiment_1002(obj,niter)
		
			v_x = -10:0.01:10;
			v_y = v_x.^2;
			plot(v_x,v_y)
			
			F = GFigure.captureCurrentFigure();	
		end	
				
		% This experiment is similar to 1002, but it illustrates how to
		% create more than one figure.
		function F = experiment_1003(obj,niter)
		
			figure(1)
			v_x = -10:0.01:10;
			v_y = v_x.^2;
			plot(v_x,v_y)			
			F(1) = GFigure.captureCurrentFigure();
			
			figure(2)
			v_x = -5:0.01:5;
			v_y = v_x.^3;
			plot(v_x,v_y)			
			F(2) = GFigure.captureCurrentFigure();
			
		end
				
		% Example of experiment that displays a figure. An GFigure is
		% constructed and returned. GSim will use this GFigure to display
		% a figure. 
		function F = experiment_1004(obj,niter)

			v_x = -10:0.01:10;
			v_y = v_x.^3;			
			F = GFigure('m_X',v_x,'m_Y',v_y);
			
		end

			
		
	end
	
	
end
