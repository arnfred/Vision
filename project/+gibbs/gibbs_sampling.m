function u = gibbs_sampling(mrf, dims, max_iter)

	% Calculate B
	B		= gibbs.get_B(mrf);

	% Initialize Burn-in
	burn	= 1;

	% Now loop through the gibbs sampling
	for i = 1:max_iter
		
		% First get z
		z		= gibbs.sample_z(mrf, u_old, dims);

		% Then get u
		u		= gibbs.sample_u(mrf, z, B, sigma, dims)

		% If we are in the burn in phase, check if we should move on
		if (burn == 1)	burn = gibbs.burn_in_p(u);

		% If not, then check what???? (Then we start denoising the image)
		else			break; end

	end
end




