
function [] = create_plots()

	% Create folder within images with a timestamp
	folder = [datestr(now, 'yyyy-mm-dd_HHMM'), '/']; 
	mkdir(['img/', folder]);

	% Load Normal data to store an overview
	[psnr u res iter burn] = load_data_sigma('mu',1,0);

	% Create normal overview:
	plot.overview(psnr, iter, burn, [folder 'normal_overview'], 'auto');

	% Create zoom of normal overview
	plot.cg_run(psnr, iter, res, 1, 0.05, [folder 'normal_zoom'], 'auto');

	% Create plot of scales with fixed sigma
end
