% Function for sampling z (a vector of mixture scales), given u (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function u = sample_u(mrf, z, B, sigma)

	% Number of rows in diagonal
	N 			= prod(mrf.imdims) * mrf.nfilters; % size of filter (82*82) times number of filters (8)

	% Calculate right hand side
	rhs			= B' * matDiag(sqrt(z)) * randn(N, 1);

	% Now calculate u
	[u, k, res]	= cg(B, z, rhs, mrf.epsilon);

	conjugate_gradient_iterations = k

end
