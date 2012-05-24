% Denoise function using gibbs sampling

function [img_denoised] = denoise()

	% Set variables
	sigma 					= 10;
	mrf 					= learned_models.cvpr_3x3_foe;
	border					= 9; % Specified in darmstadt paper
	max_iter				= 400;
	max_burnin				= 100;

	% Get images
	[img_clean img_noisy]	= denoise_init_img(sigma);

	% Pad noisy image and note image dimensions
	N_padded				= pml.support.mirror_boundary(img_noisy, border);
	u						= N_padded(:);
  	mrf.imdims				= size(N_padded);

	% Initialize statistics
	psnrs					= zeros(1, max_iters);
	ssims					= zeros(1, max_iters); % NOTE: comparison index, I figure
	mapd					= zeros(1, max_iters);
	cpu_time				= zeros(1, max_iters);

	% Get filter matrix
	B						= gibbs.get_B(mrf);

	% Run the gibbs sampling for burn i
	u						= gibbs.sample(mrf, u, B, sigma, max_burnin); % 100 is max_iter for burn_in

	for i = 1:max_iter

		% get a denoised image
		u						= img_denoise(u, B);

		% Check if the denoised image is below a certain treshold
		ssims					= stats.ssims(u, img_clean);

		% TODO: find out exactly what measurements are needed here
	end






end
