% Function for sampling z (a vector of mixture scales), given x (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function z = sample_z(mrf, x)

	% Sample z (sparse) from mrf
	z_cell		= mrf.sample_z(x(:));

	% format diagonal matrix (Cut off the end. This should be changed in mrf.sample_z)
	z			= cell2mat(z_cell);

end
