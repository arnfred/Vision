% Denoise function using gibbs sampling

function [u u_mu iter iter_mu burn_iter psnr psnr_mu] = denoise(sigma, scaling, remove_scales, max_burn, max_iter)

	% Set variables
	if (nargin < 1)	sigma				= 10; end
	if (nargin < 2)	scaling				= 1; end
	if (nargin < 3)	remove_scales		= 0; end
	if (nargin < 4)	max_burn			= 5000; end  % These are cg iterations (about 200 per gibbs iteration)
	if (nargin < 5)	max_iter			= 20000; end % Also cg iterations

	% Get images
	[img_clean img_noisy]				= denoise_init_img(sigma);
	border								= 9; % Specified in darmstadt paper
	N_padded							= pml.support.mirror_boundary(img_noisy, border);
	u									= N_padded(:);

	% set up MRF
	mrf									= learned_models.cvpr_3x3_foe;
  	mrf.imdims							= size(N_padded);
	mrf.conv_method 					= 'circular';

	% Modify the scales of the MRF
	mrf									= modify_scales(mrf, scaling, remove_scales);

	% Initialize function for getting psnr
	psnr_fun							= @(noisy) get_psnr(mrf, img_clean, noisy, border);

	% Get filter matrix
	B									= gibbs.get_B(mrf);

	% Run the gibbs for a while to burn in
	[u_b u_mu_b i_b i_mu_b p_b p_mu_b]	= gibbs.sample(mrf, u, B, sigma, max_burn, psnr_fun);

	% Now run the sampling to collect the denoised image
	[u_s u_mu_s i_s i_mu_s p_s p_mu_s]	= gibbs.sample(mrf, u_b(:,end), B, sigma, max_iter, psnr_fun);

	% Collect data for output
	u									= [u_b u_s]; 				% each iteration of the image
	u_mu								= [u_mu_b u_mu_s];			% each iteration of the mean of the image
	iter								= [i_b i_s];				% How many iterations for each image
	iter_mu								= [i_mu_b i_mu_s];			% Iterations for calculating the mean
	burn_iter							= [sum(i_b) sum(i_mu_b)];	% The total sum of burn iterations
	psnr								= [p_b p_s];				% The signal to noise ration per iteration
	psnr_mu								= [p_mu_b p_mu_s];			% The signal to noise ration per iteration for the rao blackwellisation

end
