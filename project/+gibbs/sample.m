function u = sample(mrf, dims, start_img, max_burn, sigma)

	% Default values
	if (nargin < 3)
		start_img	= randn(dims); % ~N(0,1)
		max_burn	= 100;
		sigma		= 10;
	end

	% Add path
	addpath('mat');
	addpath('inf');

	% Calculate B
	B		= gibbs.get_B(mrf, dims);

	% Initialize Burn-in and treshold
	burn	= 1, burn_treshold = 0.1;

	% Since I don't have a noise image for the pure gibbs sampling, I'll start out with a gaussian vector 
	u		= reshape(imfilter(start_img, fspecial('gaussian', 3, sigma)), prod(dims), 1);

	% Now loop through the gibbs sampling
	for i = 1:max_burn

		i
		
		% First get z
		z			= gibbs.sample_z(mrf, u);

		% Then get u
		u			= gibbs.sample_u(mrf, z, B, sigma);

		% Give me something to look at while I wait
		u_norm 		= (u - min(u(:))) / (max(u) - min(u));
		imshow(reshape(u_norm,82,82))
		drawnow

		% Calculate the energy of u
		energy(i)	= mrf.energy(u);

		% If we are in the burn in phase, check if we should move on
		if (gibbs.burn_in_p(mrf, energy, burn_treshold) == 0) break; end

	end
end




