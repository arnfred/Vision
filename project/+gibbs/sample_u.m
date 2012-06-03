% Function for sampling z (a vector of mixture scales), given u (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function [u u_mean] = sample_u(mrf, z, start_sample, B, sigma)

	% Number of rows in diagonal
	N 				= prod(mrf.imdims) * mrf.nfilters; % size of filter (82*82) times number of filters (8)

	% Calculate right hand side
	rhs				= B' * matDiag(sqrt(z)) * randn(N, 1);

	% Now calculate the mean and delta of u
	[u_mean i r]	= cg(B, z, start_sample / (sigma^2), mrf.epsilon + (1/ sigma^2));
	[u_delta k res]	= cg(B, z, rhs, mrf.epsilon + (1/ sigma^2));

	% U is equal to the mean plus the delta
	u				= u_mean + u_delta;

	conjugate_gradient_iterations_for_mean = i
	conjugate_gradient_iterations_for_delta = k

end
