% Function for sampling z (a vector of mixture scales), given x (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function z = sample_x(mrf, x, dims)

	% Number of rows in diagonal (Ideally I'll change the sampling code to sample the right amount of z's)
	N 			= prod(dims) * mrf.nfilters; % size of filter (80*80) times number of filters (8)

	% Sample z (sparse) from mrf
	z_sparse	= mrf.sample_z(x(:));

	% format diagonal matrix (Cut off the end. This should be changed in mrf.sample_z)
	z			= slice(full(cell2mat(z_sparse),1,N);

end
