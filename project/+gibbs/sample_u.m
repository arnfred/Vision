% Function for sampling z (a vector of mixture scales), given u (the image
% prior). So far this is a pretty simple one-liner, but I might later have to
% change it around a bit
function [u u_mean i i_rb psnr psnr_rb] = sample_u(mrf, z, start_sample, B, sigma, psnr_fun, psnr_fun_rb)

	% Number of rows in diagonal
	N 						= prod(mrf.imdims) * mrf.nfilters; % size of filter (82*82) times number of filters (8)

	% Calculate right hand side
	rhs						= B' * matDiag(sqrt(z)) * randn(N, 1);

	% Now calculate the mean of u
	[u_mean i_rb r psnr_rb]	= cg_prec(B, z, start_sample / (sigma^2), mrf.epsilon + (1/ sigma^2), psnr_fun_rb);

	% Prepare psnr_fun for u_delta
	psnr_fun_delta			= @(u_delta) psnr_fun(u_mean + u_delta);

	% Calculate the delta of u
	[u_delta i res psnr]	= cg_prec(B, z, rhs, mrf.epsilon + (1/ sigma^2), psnr_fun_delta);

	% U is equal to the mean plus the delta
	u							= u_mean + u_delta;

	% output
	cg_iter_for_mean 			= i_rb
	cg_iter_for_delta 			= i

end
