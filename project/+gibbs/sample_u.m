% Function for sampling z (a vector of mixture scales), given u (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function [u u_mu i i_mu psnr psnr_mu r r_mu] = sample_u(mrf, z, start_sample, B, sigma, psnr_fun, psnr_fun_mu, tolerance, use_prec)

	% Number of rows in diagonal
	N							= prod(mrf.imdims) * mrf.nfilters; % size of filter (82*82) times number of filters (8)

	% Calculate right hand side
	rhs							= B' * matDiag(sqrt(z)) * randn(N, 1);

	% Check if we want to use the preconditioner
	if (use_prec)
		% Now calculate the mean of u
		[u_mu i_mu r_mu psnr_mu]	= cg_prec(B, z, start_sample / (sigma^2), mrf.epsilon + (1/ sigma^2), psnr_fun_mu, tolerance);

		% Prepare psnr_fun for u_delta
		psnr_fun_delta				= @(u_cg) psnr_fun(u_mu + u_cg);

		% Calculate the delta of u
		[u_delta i r psnr]			= cg_prec(B, z, rhs, mrf.epsilon + (1/ sigma^2), psnr_fun_delta, tolerance);

	else % Don't use the preconditioner

		% Now calculate the mean of u
		[u_mu i_mu r_mu psnr_mu]	= cg(B, z, start_sample / (sigma^2), mrf.epsilon + (1/ sigma^2), psnr_fun_mu, tolerance);

		% Prepare psnr_fun for u_delta
		psnr_fun_delta				= @(u_cg) psnr_fun(u_mu + u_cg);

		% Calculate the delta of u
		[u_delta i r psnr]			= cg(B, z, rhs, mrf.epsilon + (1/ sigma^2), psnr_fun_delta, tolerance);
	end


	% U is equal to the mean plus the delta
	u							= u_mu + u_delta;
	psnr(end)
	psnr_mu(end)

	% output
	cg_iter_for_mean 			= i_mu
	cg_iter_for_delta 			= i

end
