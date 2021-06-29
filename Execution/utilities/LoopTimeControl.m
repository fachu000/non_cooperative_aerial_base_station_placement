classdef LoopTimeControl < handle;
	%Loop Time Control. This class allows to control the time for very
	%time-consuming loops.
	% EXAMPLE:
	% % I want to compute the cholesky decomposition of a random
	% % matrix 50 times.
	% ltc = LoopTimeControl(500);
	% for i = 1:500
	%     A = rand(1000);
	%     B = inv(A);
	%     ltc.go(i);
	% end
	
	properties
		tBegin
		tLast
		elapsedTime
		eta
		myString = [];
		totalIter;
	end
	
	methods
		function obj = LoopTimeControl(niter)
			obj.tBegin = tic;
			obj.tLast = tic;
			obj.totalIter = niter;
		end
		
		function go(obj, i)
			if toc(obj.tLast) > 0.1 || i==obj.totalIter
				obj.tLast = tic;
				obj.elapsedTime = toc(obj.tBegin);
				obj.eta = (obj.totalIter - i) * obj.elapsedTime ./ i;
				fprintf(repmat('\b', 1, length(obj.myString)));
				obj.myString = sprintf(...
					'Iteration %0.4g out of %d. ETA = %s', ...
					i, obj.totalIter, print_time(obj.eta));
				fprintf(obj.myString);
				if i==obj.totalIter
					fprintf(';; Done!\n');
				end
			end
		end
	end
	
end

