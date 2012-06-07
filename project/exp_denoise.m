function [] = exp_denoise(scaling, remove_scales, tolerance)

	% set fixed variables
	max_burn		= 3000;
	max_iter		= 7000;
	sigma			= [0.05 0.1 0.15];

	% Now for three different values of sigma
	for i = 1:3
		[u u_mean iter iter_mean burn_iter psnr psnr_rb r r_mu] = denoise(sigma(i), scaling, remove_scales, max_burn, max_iter, tolerance);
		fname = ['data/exp_sigma_',num2str(sigma(i)),'_scaling_',num2str(scaling),'_rem_scales_',num2str(remove_scales),'_maxBurn_',num2str(max_burn), '_maxIter_', num2str(max_iter), '_tolerance_', num2str(tolerance), '.mat'];
		save(fname,'u','u_mean','iter','iter_mean','burn_iter','psnr','psnr_rb','r','r_mu')
	end

	exit;
end
