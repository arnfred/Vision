function [f] = rem_scales(psnr, iter, burn_iter, name, auto)
	
	% Set some sane defaults
	if (nargin < 5) auto = 0; end
	if (nargin < 4) name = 0; end
	
	% Use multi_plot
	plot.multi_plot(psnr, iter, burn_iter, [0 2 4 6], 'Scales', [1 5000 24 26], name, auto);
end
