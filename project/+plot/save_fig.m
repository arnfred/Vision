
function [] = save_fig(f, name)
	% tightInset = get(gca, 'TightInset');
	% position(1) = tightInset(1);
	% position(2) = tightInset(2);
	% position(3) = 1 - tightInset(3)% - tightInset(1);
	% position(4) = 1 - tightInset(4)% - tightInset(2);
	% set(gca, 'Position', position);
	saveas(f, ['img/',name, '.pdf']);
end
