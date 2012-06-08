
function [] = create_plots()

	% Create folder within images with a timestamp
	folder = [datestr(now, 'yyyy-mm-dd_HHMM'), '/']; 
	mkdir(['img/', folder]);

	% Set larger font-size
	set(0,'DefaultAxesFontSize',16)

	% Load Normal data with sigma variations to store an overview with preconditioner
	[psnr1 u1 res1 iter1 burn1] = load_data_sigma('mu',1,0);

	% Create normal overview:
	plot.sigma(psnr1, iter1, burn1, [folder 'sigma_with_prec'], 'auto');

	% Load Normal data with sigma variations to store an overview without preconditioner
	[psnr1_1 u1_1 res1_1 iter1_1 burn1_1] = load_data_sigma('mu',1,0,0.001,0);

	% Create normal overview:
	plot.sigma(psnr1_1, iter1_1, burn1_1, [folder 'sigma_without_prec'], 'auto');

	% Create zoom of normal overview
	plot.cg_run(psnr1, iter1, res1, 1, 0.05, [folder 'normal_zoom'], 'auto');

	% Export some images to use
	plot.get_picture_at_psnr(psnr1, iter1, u1, 30, [folder 'img_30_psnr'], 'auto');
	plot.get_picture_at_psnr(psnr1, iter1, u1, 28, [folder 'img_28_psnr'], 'auto');
	plot.get_picture_at_psnr(psnr1, iter1, u1, 26, [folder 'img_26_psnr'], 'auto');
	plot.get_picture_at_psnr(psnr1, iter1, u1, 24, [folder 'img_24_psnr'], 'auto');
	plot.get_picture_at_psnr(psnr1, iter1, u1, 22, [folder 'img_22_psnr'], 'auto');




	% Load Normal data with rem_scales variation
	[psnr2 u2 res2 iter2 burn2] = load_data_rem_scales('mu',1,0.1);
	
	% Create plot of scales with fixed sigma
	plot.rem_scales(psnr2, iter2, burn2, [folder 'rem_scales_overview'], 'auto');





	% Load data for scale variation
	[psnr3 u3 res3 iter3 burn3] = load_data_scaling('mu',0,0.1,[1 0.8 0.6 0.4]);
	% Plot the same data we just loaded
	plot.multi_plot(psnr3, iter3, burn3, [1 0.8 0.6 0.4], 'Scale Factor', [1 5000 23 26], (1:4), [folder 'scaling_overview'], 'auto');
	% [psnr3 u3 res3 iter3 burn3] = load_data_scaling('mu',0,0.1,[1 0.95 0.9 0.85]);
	% Plot the same data we just loaded
	% plot.multi_plot(psnr3, iter3, burn3, [1 0.95 0.9 0.85], 'Scale Factor', [1 5000 23 26], (1:4), [folder 'scaling_top'], 'auto');


	% Get data and make plot for tolerance
	[psnr4 u4 res4 iter4 burn4] = load_data_tolerance('mu',0.1);
	plot.multi_plot(psnr4, iter4, burn4, [0.0001 0.001 0.01 0.1], 'Tolerance', [1 5000 24 26], fliplr(1:4), [folder 'tolerance_overview'], 'auto');
	plot.multi_plot(psnr4, iter4, burn4, [0.0001 0.001 0.01 0.1], 'Tolerance', [1 5000 24 26], (1:4), [folder 'tolerance_flipped'], 'auto');

	% Get data and make plot for tolerance
	[psnr6 u6 res6 iter6 burn6] = load_data_tolerance('mu',0.05);
	plot.multi_plot(psnr6, iter6, burn6, [0.0001 0.001 0.01 0.1], 'Tolerance', [1 5000 28 30.5], fliplr(1:4), [folder 'tolerance_low_sigma'], 'auto');

	% Present the full scaling data (this takes a while to load)
	[psnr5 u5 res5 iter5 burn5] = load_data_scaling('mu',0,0.1);
	plot.heatmap(psnr5, iter5, [folder 'Scaling_heatmap'], 'auto');



end
