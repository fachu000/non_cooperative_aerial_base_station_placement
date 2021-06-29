
		
		function Y = project_2D_along( w , X )
			%
			% Projection of a collection of points in the 2D plane 
			% 
			% INPUT:
			%  X     : is a p x N matrix whose columns are the points to be 
			%         projected on R^2
			%  w     : is a p x 1 vector which determines one of the directions
			%         of the projection. The other direction is that orthogonal
			%         to "w" that maximizes the variability of the columns of
			%         "X". 
			% OUTPUT:
			%  Y     : is a 2 x N vector with the coordinates of the projections
			%
		    % Every column of "Y", say "y" is computed as:
			%
			%             y = [ w' ; v']*x
			%
			% Where "x" is the solution of the problem
			%
			%          maximize     || v' * X ||_2
			%          s.t.         w'*v = 0
			%                       || v || <= 1
			%
			% "w" is normalized to unit norm if the input does not satisfy
			% this property.
			
			p = size(w);
			w = w/norm(w);
			N = size(X,2);
			
			assert(size(X,1)==p);
			
			
			% the following can also be computed directly in terms of eigvals
			cvx_begin 
			    variable v(p)    
			    maximize( norm(v'*X) )
				such that
				   w'*v == 0;
				   norm(v)<=1;
			cvx_end
						             
			Y  = [ w' ; v']*X;
			
		end
		