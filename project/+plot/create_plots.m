
function [] = create_plots()

	% Create folder within images with a timestamp
	folder = [datestr(now, 'yyyy-mm-dd_HHMM'), '/']; 
	mkdir(['img/', folder]);


	% Load Normal data with sigma variations to store an overview
	[psnr1 u1 res1 iter1 burn1] = load_data_sigma('mu',1,0);

	% Create normal overview:
	plot.sigma(psnr1, iter1, burn1, [folder 'sigma_overview'], 'auto');

	% Create zoom of normal overview
	plot.cg_run(psnr1, iter1, res1, 1, 0.05, [folder 'normal_zoom'], 'auto');




	% Load Normal data with rem_scales variation
	[psnr2 u2 res2 iter2 burn2] = load_data_rem_scales('mu',1,0.1);
	
	% Create plot of scales with fixed sigma
	plot.rem_scales(psnr2, iter2, burn2, [folder 'rem_scales_overview'], 'auto');


	% Load data for scale variation
	[psnr3 u3 res3 iter3 burn3] = load_data_scaling('mu',0,0.1,[1 0.8 0.6 0.4]);
	% Plot the same data we just loaded
	plot.multi_plot(psnr3, iter3, burn3, [1 0.8 0.6 0.4], 'Scale Factor', [1 5000 23 26], [folder 'scaling_overview'], 'auto');

	% What more before I can go home and sleep a little?

	% - Make plot of cg cut-off (or at least make it ready for tomorrow)
	% - Think about how to present the full scaling data
end
