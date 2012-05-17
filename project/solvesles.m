
% Matlab script for solving a bunch of linear systems based on the fields of
% experts mrf % This is very expensive when we use valid fft and not circ
function [result elapsed res] = solvesles()

  % Dimensions
  dims = [82 82];

  % Add path
  addpath('mat');
  addpath('inf');

  % Get mrf
  mrf = learned_models.cvpr_3x3_foe;
  nfilters = mrf.nfilters

  % Start out with an empty B matrix. NOTE: mat/@mat/diagFAtAFt.m was commented
  % out at line 16-17 to prevent it from returning and empty matrix
  B = mat([0 prod(dims)], 'circ');

  % For each filter, create a fake convolution matrix
  for i=1:mrf.nfilters
      Bi = matConv2(mrf.filter(i), dims, 'circ');
      B = vertcat(B,Bi); % B is (51200 + 6724) x 6724
  end

  % Load Z matrices
  load('z.mat');

  % Now we should have a cell array called zes in scope
  [~, number_of_zs] = size(zes);

  % Set up return matrix as the last 4 matrices:
  result = zeros(dims(1), dims(2), 3);

  % Number of rows in diagonal
  N = prod(dims) * nfilters; % size of filter (80*80) times number of filters (8)

  % initialize elapsed
  for i = 1:3 elapsed{i} = 0; end

  % Iterate over all Z's, solving the systems one by one
  for i = 1:number_of_zs

	  % Let's do this the right way, this time

	  % Get diagonal matrix
	  diagonal		= slice(full(cell2mat(zes(i))),1,N);

	  % Calculate right hand side
	  rhs			= B' * matDiag(sqrt(diagonal)) * randn(N, 1);

	  % Now calculate u
	  [u, k, res]	= cg(B, diagonal, rhs, mrf.epsilon);

  end

  result(:,:,1)	= reshape(u{1}, dims);
  result(:,:,2)	= reshape(u{2}, dims);
  result(:,:,3)	= reshape(u{3}, dims);

	  
end
