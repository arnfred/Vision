
function [] = save_fig(f, name)
	addpath('export_fig');
	full_path = ['img/',name, '.pdf'];
	% tightInset = get(gca, 'TightInset');
	% position(1) = tightInset(1);
	% position(2) = tightInset(2);
	% position(3) = 1 - tightInset(3)% - tightInset(1);
	% position(4) = 1 - tightInset(4)% - tightInset(2);
	% set(gca, 'Position', position);

	
	figureHandle = gcf;
	%# make all text in the figure to size 14 and bold
	set(findall(figureHandle,'type','text'),'fontSize',16)
    
    saveas(f, full_path);
	%export_fig(full_path);
end
