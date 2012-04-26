
% Inputs:
%	A, b	The system to be solved: Ax = b. I expect A to be a Mat from Hannes package

% Outputs
%	x		The resulting vector
%	k		The amount of iterations done
%	res		The distance in the end

% Adapted from "An Introduction to the Conjugate Gradient Method Without the
% Agonizing Pain (Jonathan Richard Shewchuk)".

function [x, i, delta_n, res] = conjugate_gradients(A, b, Prec)

	% Get sizes
	[n, m]	= size(A);

	% Initialize x. Maybe it should be some other value?
	x		= zeros(m,1);

	% Initialize i (iterations), max_i (maximum amount of iterations), epsilon
	% (error treshold)
	i		= 0;
	max_i	= 499;
	epsilon	= 1e-3;

	% Initialize Residual It should be b - Ax, but since x is initialized to 0,
	% I'm skipping that. If it's changed later, that should be fixed though.
	r		= b; % - Ax;

	% Initialize direction. The preconditioner better be easy to invert, or we
	% are wasting our time.
	d		= Prec(r);

	% Initialize deltas
	delta_n	= r'*d;
	delta_o	= delta_n;

	% Vector holding the Frobenius norm of r
	res		= zeros(max_i+1,1);
	res(1)	= norm(r,'fro');
	nb 		= norm(b,'fro');                        % |b|

	% Loop through the thing
	for i=1:max_i

		% Calculate Ad to use for later
		Ad		= A*d;
		
		% Calculate the amount we are adding of the new direction
		alpha	= delta_n/(d'*Ad);

		% Our next guess for x
		x		= x + alpha*d;

		% Should the residual be stabilized?
		if (mod(i,50) == 0) r = b - A*x;
		else				r = r - alpha*Ad;
		end
		
		% Update the norm vector
		res(i+1) = norm(r,'fro');

		% Calculate M^-1 * r to use for later
		invMr	= Prec(r);
		
		% Update deltas
		delta_o	= delta_n;
		delta_n	= r'*invMr;

		% Update beta
		beta	= delta_n / delta_o;

		% Update to the next direction
		d		= invMr + beta*d;

		% Update index and quit if we are within treshold
		if (delta_n < epsilon^2*delta_o) break; end
		if (res(i+1) < epsilon*nb); break; end
	end




