function u = gibbs_sampling(mrf, u_old, B, dims)

	% First get z
	z		= gibbs.sample_z(mrf, u_old, dims);

	% Then get u (well this is simple...)
	u		= gibbs.sample_u(mrf, z, B, sigma, dims)

