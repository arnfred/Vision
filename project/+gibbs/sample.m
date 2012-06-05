function [samples samples_mu iterations_mu iterations psnrs psnrs_mu] = sample(mrf, start_sample, B, sigma, max_iter, psnr_fun)

	% Initialize statistics
	psnrs					= []; %zeros(1, max_iter);
	psnrs_mu				= []; %zeros(1, max_iter);
	iterations				= []; %zeros(1, max_iter);
	iterations_mu			= []; %zeros(1, max_iter);
	samples					= []; %zeros(numel(start_sample), max_iter);
	samples_mu				= []; %zeros(numel(start_sample), max_iter);

	% Since I don't have a noise image for the pure gibbs sampling, I'll start out with a gaussian vector 
	% u		= reshape(imfilter(start_sample, fspecial('gaussian', 3, sigma)), prod(dims), 1);
	u				= start_sample;
	u_denoised		= start_sample;
	u_denoised_mu	= start_sample;

	% To make the image appear
	f 				= figure;
	k				= 1;

	% Now loop through the gibbs sampling
	while sum(iterations) < max_iter

		k
		
		% First get z (In the darmstadt paper they use the initial sample)
		z						= gibbs.sample_z(mrf, u);

		% setup psnr for conjugate gradients
		psnr_cg					= @(u) psnr_fun((u + (i-1) * u_denoised) / 1);
		psnr_cg_mu				= @(u) psnr_fun((u_mu + (i-1) * u_denoised_mu) / 1);

		% Then get u and u_mean
		[u u_mu i i_mu p p_mu]	= gibbs.sample_u(mrf, z, start_sample, B, sigma, psnr_cg, psnr_cg_rb);

		% Give me something to look at while I wait
		imagesc(reshape(u,82,82))
		colormap('gray')
		drawnow
		% pause

		% Calculate the denoised image as a moving average
		u_denoised				= (u + (k-1)*u_denoised) / k;
		u_denoised_mu			= (u_mu + (k-1)*u_denoised_mu) / k;

		% Calculate the psnr and store it and the iterations
		% psnr					= psnr_fun(u_denoised)
		% psnr_rb					= psnr_fun(u_denoised_rb)
		psnrs					= [psnrs p];
		psnrs_mu				= [psnrs_mu p_mu];
		iterations				= [iterations i];
		iterations_mu			= [iterations_mu i_mu];
		samples					= [samples u];
		samples_mu				= [samples_mu u_mu];

		% For following along
		i
		sum_iter				= sum(iterations) + sum(iterations_mu)

		% Update i
		k						= k + 1;

	end

end




