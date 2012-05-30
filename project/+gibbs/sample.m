function u = sample(mrf, start_sample, B, sigma, max_iter)

	% Add path
	addpath('mat');
	addpath('inf');

	% Initialize Burn-in and treshold
	burn_treshold	= 0.05;
	max_burn		= 100;

	% Since I don't have a noise image for the pure gibbs sampling, I'll start out with a gaussian vector 
	% u		= reshape(imfilter(start_sample, fspecial('gaussian', 3, sigma)), prod(dims), 1);
	u				= start_sample;

	% Now loop through the gibbs sampling
	for i = 1:max_iter

		i
		
		% First get z
		z				= gibbs.sample_z(mrf, u);
		size(z)
		plot(sort(z))
		pause

		% Then get u
		u				= gibbs.sample_u(mrf, z, B, sigma);

		% Give me something to look at while I wait
		imagesc(reshape(u,82,82))
		colormap('gray')
		drawnow
		pause

		% Calculate the energy of u
		energy(i)		= mrf.energy(u);

		% If we are in the burn in phase, check if we should move on
		if (gibbs.burn_in_p(mrf, energy, burn_treshold) == 0) break; end

	end

	% TODO: return average of the last half of the samples
end




