function [f] = overview(psnr, iter, burn_iter, name)

	% Get new figure
	f = figure;
	hold on

	% Find the maximum height of the figure
	padding		= (max(psnr) - min(psnr)) / 5;
	max_h		= max(psnr) + padding;
	min_h		= min(psnr) - padding;
	iter_sum	= cumsum(iter);

	% For each second iter, paint a gray box
	for i = 1:2:(numel(iter) - 1)
		ha = area([iter_sum(i) iter_sum(i+1)], [max_h max_h]);
		set(ha(1),'FaceColor',[.95 .98 .98],'EdgeColor',[.95 .98 .98])
	end

	% Now plot the psnr
	plot(psnr,'Color',[.1 .2 .7])
	axis([1 iter_sum(end) min_h max_h])
	hold off

	% Plot the burn-in cutoff
	line([burn_iter burn_iter], [min_h max_h], 'Color', [0.7 0.05 0.05]);

	% Set background to white
	set(gcf, 'Color', 'w');

	% Make axis stay on top of areas
	set(gca,'Layer','top')

	% If name is set, export
	if (exist('name') == 1) 
		plot.save_fig(f, name); 
	end

	% Add names and legend
	xlabel('Conjugate Gradients Iterations');
	ylabel('PSNR');

	% pause and close
	pause();
	close();
end
