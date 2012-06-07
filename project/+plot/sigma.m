function [f] = sigma(psnr, iter, burn_iter, name, auto)

	% Set some sane defaults
	if (nargin < 5) auto = 0; end
	if (nargin < 4) name = 0; end
	
	% Use multi_plot
	plot.multi_plot(psnr, iter, burn_iter, [0.05 0.1 0.15], 'Scales', [1 5000 16 32], [1 2 3], name, auto);

end
