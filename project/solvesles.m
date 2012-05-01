
% Matlab script for solving a bunch of linear systems based on the fields of
% experts mrf
function [result elapsed res] = solvesles()

  % Dimensions
  dims = [82 82];

  % Add path
  addpath('mat');
  addpath('inf');

  % Get mrf
  mrf = learned_models.cvpr_3x3_foe;
  nfilters = mrf.nfilters

  % FOR TESTING: TODO: change back
  % nfilters = 1;

  % Express B in terms of the filters of the mrf
  % B = ones(0,prod(dims));

  % Start with identity matrix
  B = matDiag(ones(prod(dims),1));

  % For each filter, create a fake convolution matrix
  for i=1:mrf.nfilters
      Bi = matConv2(mrf.filter(i), dims, 'valid');
      B = vertcat(B,Bi); % B is 51200 x 6724
  end

  % FOR TESTING: TODO: change back
  % B = matConv2(mrf.filter(1), dims, 'circ');

  % Load Z matrices
  load('z.mat');

  % Now we should have a cell array called zes in scope
  [~, number_of_zs] = size(zes);

  % Set up return matrix as the last 4 matrices:
  result = zeros(dims(1), dims(2), 3);

  % Number of rows in diagonal (dims - 2 because the 'valid' fft transform gives back a smaller matrix) 
  % N = prod(dims - 2) * nfilters; % size of filter (80*80) times number of filters (8)
  N = prod(dims - 2) * nfilters + prod(dims); % Add another 6724 for the diagonal
  % TODO: this is only for testing
  % N = prod(dims) * nfilters; % size of filter (80*80) times number of filters (8)

  % initialie elapsed
  for i = 1:3 elapsed{i} = 0; end

  % Iterate over all Z's, solving the systems one by one
  for i = 1:number_of_zs
  % for i = 1:5
	  
	  % Get diagonal matrix
	  origDiagonal	= full(cell2mat(zes(i)));

	  % Cut off identity matrix saved at the end
	  % diagonal		= matDiag(origDiagonal(1:N));
	  diagonal		= matDiag(origDiagonal);
	  % sqrtDiagonal	= matDiag(sqrt(origDiagonal(1:N)));
	  sqrtDiagonal	= matDiag(sqrt(origDiagonal));
     
	  % Calculate right hand side
	  rhs			= B' * sqrtDiagonal * randn(N, 1); % n0 in notes

	  % Prepare preconditioners
	  A			= B'*diagonal*B;
	  diag_A		= diag(A);
	  d 			= mean(diag(diagonal)) * diagFAtAFt(B,'cheap');
	  F 			= matFFTNmask(true(dims)); % DFT matrix
	  M{1}			= @(r) r;
	  M{2}			= @(r) r ./ diag_A;
	  M{3}			= @(r) [F']*((F*r)./d(:));

	  % Loop over the different preconditioning schemes
	  for (j = 1:3) 

		  tstart{j}					= cputime;
		  [u{j} iter{j} d res{j}]	= conjugate_gradients(A, rhs, M{j});
		  elapsed{j}				= elapsed{ j } + cputime - tstart{ j };

	  end

  end

  result(:,:,1)	= reshape(u{1}, dims);
  result(:,:,2)	= reshape(u{2}, dims);
  result(:,:,3)	= reshape(u{3}, dims);

	  
end
