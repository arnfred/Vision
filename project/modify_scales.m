% Modifies the scales of an MRF

function mrf = modify_scales(mrf, scaling, remove_scales)

if (mod(remove_scales,2) == 1) error('The amount of scales to remove must be divisible with two'); end

	% For each expert, adjust the scales and drop the highest ones
	nexperts = numel(mrf.experts);
	for i = 1:nexperts
		experts					= mrf.experts{i};
		scales					= [-9, -7 -5:5, 7, 9];
		remove					= remove_scales/2;

		% Scale
		scales					= scales * scaling;

		% Remove
		scales					= scales((remove+1):end-remove);

		% Get new scales
		new_scales				= exp(scales);

		% Update variables
		mrf.experts{i}.scales	= exp(scales);
		mrf.experts{i}.weights	= experts.weights((remove+1):end-remove);
		% mrf.experts{i}.nscales	= experts.nscales - remove_scales;
	end

end
