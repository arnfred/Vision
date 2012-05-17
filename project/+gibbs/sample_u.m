% Function for sampling z (a vector of mixture scales), given u (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function u = sample_u(mrf, z, B, sigma, dims)

	% Number of rows in diagonal
	N 			= prod(dims) * mrf.nfilters; % size of filter (80*80) times number of filters (8)

	% Get diagonal matrix
	diagonal	= slice(full(cell2mat(z),1,N);

	% Calculate right hand side
	rhs			= B' * matDiag(sqrt(z)) * randn(N, 1);

	% Now calculate u
	[u, k, res]	= cg(B, z, rhs, mrf.epsilon);

end
