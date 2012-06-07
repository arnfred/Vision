
function [f] = heatmap(psnr, iter, name, auto)

	% Get some variables
	if (nargin < 4) auto = 0; end
	if (nargin < 3) name = 0; end

	elements		= numel(psnr);
	span			= 100;
	max_elem		= 5000;

	% Get min and max
	top = 0;
	bottom = 9999999;
	for i = 1:elements
		top			= max(top, max(psnr{i}));
		bottom		= min(bottom, min(psnr{i}));
	end

	% Initialize my matrix
	indices			= floor(max_elem/span);
	matrix			= ones(elements,floor(max_elem/span)) * bottom;

	% For each psnr
	for k = 1:elements

		% Get the element
		e = psnr{k};

		% For each 100 iterations of cg find the max value
		for i = 1:indices
			range = ((span * (i - 1)) + 1):span*i;

			% Check if we are within range of the element
			if (range(end) < numel(e))

				% Get the maximum point between i and i - span
				max_point	= max(e(range));

				% Add this point to the matrix
				matrix(k,i)	= max_point;

			else
				matrix(k,i) = matrix(k,i-1);
			end

		end

	end

	% Create a colormap
	color1			= [.95 .99 .99].^20;
	color2			= [.99 .95 .99].^20;
	grid1			= meshgrid(color1,1:40);
	grid2			= meshgrid(color2,1:40);
	part1			= grid1' * diag((1:40)/20);
	part2			= grid2' * diag(fliplr(1:40)/20);
	map				= part1' + part2';

	% Normalize
	map				= map / max(map(:));

	% Now create a figure
	f = figure;

	imagesc(matrix);
	colormap(map);
	cb = colorbar;

	% Set background to white
	set(gcf, 'Color', 'w');

	% Set the right axis labels
	xTicks = get(gca,'XTick')*span;  % Get the x axis ticks and multiply them with the span
	yTicks = 1 - (get(gca,'YTick') - 1) * 0.05;  % Get the x axis ticks and multiply them with the span
	set(gca,'XTickLabel',num2str(xTicks.'));
	set(gca,'YTickLabel',num2str(yTicks.'));

	% Add axis labels
	xlabel('Conjugate Gradients Iterations');
	ylabel('Scale Factor');

	% Add colorbar label
	set(get(cb,'ylabel'),'string','psnr')

	% If name is set, export
	if (name ~= 0) 
		plot.save_fig(f, name); 
	end

	% pause and close
	if (auto == 0) pause(); end
	close();

end
