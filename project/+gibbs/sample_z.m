% Function for sampling z (a vector of mixture scales), given x (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function z = sample_z(mrf, x)

	% Sample z (sparse) from mrf
	z_cell		= mrf.sample_z(x(:));

	% Scale the scales
	for i = 1:mrf.nfilters
		nexperts 			= mrf.nexperts;
		expert_precision 	= mrf.experts{min(i,nexperts)}.precision;
		z_scaled{i} 		= expert_precision * z_cell{i}(:);
	end

	% format diagonal matriu_thiersx
	z			= cell2mat(z_scaled');

end
