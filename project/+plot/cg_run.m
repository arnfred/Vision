% Plots the run of one cg iteration
function [] = cg_run(psnr, iter, res, iteration, s, name, auto)

	% Set some sane defaults
	if (nargin < 5) s = 0.1; end
	if (nargin < 4) iteration = 10; end

	% Get new figure
	f = figure;
	hold on

	% Get index
	index = 0;
	if (s == 0.05) index = 1; end
	if (s == 0.1) index = 2; end
	if (s == 0.15) index = 3; end
	if (index == 0) error(['There is no data for the sigma value: ',num2str(s)]); end

	% Get data
	iter_sum		= [1 cumsum(iter{index})];
	ind				= (iter_sum(iteration):iter_sum(iteration + 4));
	psnr			= psnr{index}(ind);
	res				= res{index}(ind);
	colors_line		= [.005 .7 .7; .7 .005 .7; .7 .7 .005; .3 .3 .3];
	colors_area		= [.95 .99 .99; .99 .95 .99; .99 .99 .95; .95 .95 .95];

	size(ind)

	% Now plot the psnr
	x = 1:numel(res);
	[ax h1 h2] = plotyy(x, psnr, x, res, 'semilogy', 'semilogy')
	set(ax,{'ycolor'},{colors_line(index,:);colors_line(4,:)})
	set(h2 , 'Color', colors_line(4,:), 'LineStyle', '-');
	set(h1 , 'Color', colors_line(index,:));

	% Add axis labels
	axes(ax(1)); xlabel('Conjugate Gradients Iterations'); ylabel(['psnr  for sigma of ', num2str(s)]);
	axes(ax(2)); ylabel('(Norm of residual) / (norm of right hand side)');

	% Add area under line
	axes(ax(1));
	c = colors_area(index,:);
	H = area(ax(1), x, psnr);
	h = get(H,'children');
	set(H,'FaceColor',c,'EdgeColor',c)
	line(x, psnr, 'Color', colors_line(index,:));

	% Set x-axis (else there are two overlapping
	linkaxes(ax, 'x')
	set(ax(1), 'XLim', [1 numel(res)], 'YLim', [16 32]);

	% Specify y-axis
	set(ax(2), 'YLim', [10e-5 10e3]);


	hold off

	% Set background to white
	set(gcf, 'Color', 'w');

	% Make axis stay on top of areas
	axes(ax(1));
	set(gca,'Layer','top')
	axes(ax(2));
	set(gca,'Layer','top')

	% If name is set, export
	if (exist('name') == 1) 
		plot.save_fig(f, name); 
	end

	% pause and close
	if (exist('auto') == 0) pause(); end
	close();
end
