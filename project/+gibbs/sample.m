function [samples samples_mu iterations iterations_mu psnr psnr_mu res res_mu burn_iter] = sample(mrf, start_sample, B, sigma, max_burn, max_iter, psnr_fun, border, tolerance, use_prec)

	% Initialize statistics
	psnr					= []; %zeros(1, max_iter);
	psnr_mu					= []; %zeros(1, max_iter);
	iterations				= []; %zeros(1, max_iter);
	iterations_mu			= []; %zeros(1, max_iter);
	samples					= []; %zeros(numel(start_sample), max_iter);
	samples_mu				= []; %zeros(numel(start_sample), max_iter);
	res						= []; %zeros(numel(start_sample), max_iter);
	res_mu					= []; %zeros(numel(start_sample), max_iter);

	% Since I don't have a noise image for the pure gibbs sampling, I'll start out with a gaussian vector 
	% u		= reshape(imfilter(start_sample, fspecial('gaussian', 3, sigma)), prod(dims), 1);
	u				= start_sample;
	u_denoised		= start_sample;
	u_denoised_mu	= start_sample;

	% To make the image appear
	% f 			= figure;
	k				= 1;

	% Now loop through the gibbs sampling
	while (sum(iterations) + sum(iterations_mu)) < max_iter

		% Check if we are done burning in which case the running average is reset
		if (sum(iterations) + sum(iterations_mu) > max_burn) 
			k 								= 1; 
			burn_iter						= [sum(iterations) sum(iterations_mu)];
			max_burn						= 9999999999999; % Something big
		end

		k
		
		% First get z (In the darmstadt paper they use the initial sample)
		z								= gibbs.sample_z(mrf, u);

		% setup psnr for conjugate gradients
		psnr_cg							= @(u_cg) psnr_fun((u_cg + (k-1) * u_denoised) / k);
		psnr_cg_mu						= @(u_cg_mu) psnr_fun((u_cg_mu + (k-1) * u_denoised_mu) / k);

		% Then get u and u_mean
		[u u_mu i i_mu p p_mu r r_mu]	= gibbs.sample_u(mrf, z, start_sample, B, sigma, psnr_cg, psnr_cg_mu, tolerance, use_prec);

		% Calculate the denoised image as a moving average
		u_denoised						= (u + (k-1)*u_denoised) / k;
		u_denoised_mu					= (u_mu + (k-1)*u_denoised_mu) / k;

		% Give me something to look at while I wait
		% imshow(get_img(u_denoised_mu, border, mrf.imdims));
		% colormap('gray');
		% drawnow

		% Calculate the psnr and store it and the iterations
		psnr							= [psnr p];
		psnr_mu							= [psnr_mu p_mu];
		iterations						= [iterations i];
		iterations_mu					= [iterations_mu i_mu];
		samples							= [samples u];
		samples_mu						= [samples_mu u_mu];
		res								= [res; r(2:end)];
		res_mu							= [res_mu; r_mu(2:end)];

		% For following along
		sum_iter						= sum(iterations) + sum(iterations_mu)

		% Update k
		k						= k + 1;

	end

end




