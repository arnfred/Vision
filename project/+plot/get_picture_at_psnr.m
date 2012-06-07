
function [] = get_picture_at_psnr(psnr, iter, u, target, name, auto)

	% Get some variables
	if (nargin < 6) auto = 0; end
	if (nargin < 5) name = 0; end

	% Vertcat the arrays
	psnr	= [horzcat(psnr{:}) 0];
	iter	= horzcat(iter{:});
	u		= horzcat(u{:});
	i_sum	= cumsum(iter);

	% Now flip through psnr until we go above the target
	for i = 1:(numel(i_sum)-1)

		k 		= i_sum(i);
		k_next	= i_sum(i+1);
		disp(['at ', num2str(k), ' we have a psnr of ', num2str(psnr(k))]);

		% Check if target is higher than sample and lower than next sample
		if ((psnr(k) < target && psnr(k_next) > target) | ((psnr(k) - 0.1) < target && (psnr(k) + 0.1) > target))

			% Get image
			img		= get_img(u(:,i), 9, [82 82]);
			f		= imagesc(img/max(img(:)));
			colormap('gray')
			set(gca,'xtick',[], 'xticklabel',{})
			set(gca,'ytick',[], 'yticklabel',{})

			% Set background to white
			set(gcf, 'Color', 'w');

			% If name is set, export
			if (name ~= 0) 
				plot.save_fig(f, name); 
			end

			% pause and close
			if (auto == 0) pause(); end
			close();

			% Exit
			return
		end
	end

	disp(['No image with psnr around ', num2str(target), ' was found']);

end
