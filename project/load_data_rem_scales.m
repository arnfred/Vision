% Function for loading all sigma values for a certain scaling and rem_scales

function [psnr_out u_out res_out iter_out burn_out] = load_data_rem_scales(type, scaling, sigma)

	% Set boolean to 1 if we want to load the rao blackwellised data
	if (type == 'rb' | type == 'mu' | type == 'mean') mu = 1;
	else mu = 0; end

	% Initialize values
	psnr_out		= {};
	u_out			= {};
	res_out			= {};
	iter_out		= {};
	burn_out		= {};

	% Sigma and iter values as taken from exp_denoise
	rem_scales		= [0 2 4 6];
	max_burn		= 3000;
	max_iter		= 7000;

	% For each sigma value we load the file and collect the appropriate values
	for i = 1:numel(rem_scales)

		% Now find a filename that fits
		fname = ['data/exp_sigma_',num2str(sigma),'_scaling_',num2str(scaling),'_rem_scales_',num2str(rem_scales(i)),'_maxBurn_',num2str(max_burn), '_maxIter_', num2str(max_iter),'.mat'];

		% Check if it exists
		if (exist(fname) ~= 2) disp(['Could not load file with name: ',fname]); break; end
		load(fname);

		% Assign appropriate values
		if (mu == 1)
			psnr_out{i}			= psnr_rb;
			u_out{i}			= u_mean;
			res_out{i}			= r_mu;
			iter_out{i}			= iter_mean;
			burn_out{i}			= burn_iter(2);
		else
			psnr_out{i}			= psnr;
			u_out{i}			= u;
			res_out{i}			= r;
			iter_out{i}			= iter;
			burn_out{i}			= burn_iter(1);
		end
	end

end
