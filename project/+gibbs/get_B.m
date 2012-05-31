
% Returns a filter matrix based on the filters
function B = get_B(mrf)

	% Add path
	addpath('mat');
	addpath('inf');

	% Start out with an empty B matrix. NOTE: mat/@mat/diagFAtAFt.m was commented
	% out at line 16-17 to prevent it from returning and empty matrix
	B = mat([0 prod(mrf.imdims)], 'circ');

	% For each filter, create a fake convolution matrix
	for i=1:mrf.nfilters
		Bi = matConv2(mrf.filter(i), mrf.imdims, 'circ');
		B = vertcat(B,Bi); % B is (51200 + 6724) x 6724
	end

end
