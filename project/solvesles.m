
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

  % FOR TESTING: TODO: change back
  % nfilters = 1;

  % Express B in terms of the filters of the mrf
  % B = ones(0,prod(dims));

  % We don't add the identity matrix here, this must be done later
  % % Start with identity matrix
  % B = matDiag(ones(prod(dims),1));
  B = mat([0 prod(dims)], 'circ');

  % For each filter, create a fake convolution matrix
  for i=1:mrf.nfilters
      Bi = matConv2(mrf.filter(i), dims, 'circ');
      B = vertcat(B,Bi); % B is (51200 + 6724) x 6724
  end

  % FOR TESTING: TODO: change back
  % B = matConv2(mrf.filter(1), dims, 'valid');

  % Load Z matrices
  load('z.mat');

  % Now we should have a cell array called zes in scope
  [~, number_of_zs] = size(zes);

  % Set up return matrix as the last 4 matrices:
  result = zeros(dims(1), dims(2), 3);

  % Number of rows in diagonal (dims - 2 because the 'valid' fft transform gives back a smaller matrix) 
  N = prod(dims - 2) * nfilters; % size of filter (80*80) times number of filters (8)
  % N = prod(dims - 2) * nfilters + prod(dims); % Add another 6724 for the diagonal
  % TODO: this is only for testing
  % N = prod(dims) * nfilters; % size of filter (80*80) times number of filters (8)

  % initialize elapsed
  for i = 1:3 elapsed{i} = 0; end

  % Iterate over all Z's, solving the systems one by one
  for i = 1:number_of_zs

	  % Let's do this the right way, this time

	  % Get diagonal matrix
	  diagonal		= slice(full(cell2mat(zes(i))),1,N);

	  % Calculate right hand side
	  rhs			= B' * matDiag(sqrt(diagonal)) * randn(N, 1);

	  % Declare a few parameters (should idealically be set earlier)
	  identity		= matDiag(ones(1,prod(dims)));
	  kmax			= 500; 					% Run for a maximum of kmax iterations
	  c				= zeros(size(rhs)); 	% The starting vector
	  tol			= 1e-3; 				% Right now this isn't happening anyway
	  prec			= 1;					% Please use preconditioning

	  % Now calculate u
	  %[u, k, res]	= linsolve_lcg(B, diagonal, identity, mrf.epsilon, rhs, kmax, c, tol, prec);
	  [u, k, res]	= cg(B, diagonal, rhs, mrf.epsilon);





  % for i = 1:5
	  

	  % Cut off identity matrix saved at the end
	  % diagonal		= matDiag(origDiagonal(1:N));
	  % sqrtDiagonal	= matDiag(sqrt(origDiagonal(1:N)));
	  %  sqrtDiagonal	= matDiag(sqrt(origDiagonal(1:N)));
     
	  %  % Calculate right hand side
	  %  rhs			= B' * sqrtDiagonal * randn(N, 1); % n0 in notes

	  %  % Prepare preconditioners
	  %  A				= B'*diagonal*B + eye(prod(dims))*mrf.epsilon;
	  %  diag_A		= diag(A);
	  %  d 			= mean(diag(diagonal)) * diagFAtAFt(B,'cheap') + mrf.epsilon; % This is very expensive when we use valid fft and not circ
	  %  F 			= matFFTNmask(true(dims)); % DFT matrix
	  %  M{1}			= @(r) r;
	  %  M{2}			= @(r) r ./ diag_A;
	  %  M{3}			= @(r) [F']*((F*r)./d(:));

	  %  % Loop over the different preconditioning schemes
	  %  for (j = 1:3) 

	  %      tstart{j}					= cputime;
	  %      [u{j} iter{j} d res{j}]	= conjugate_gradients(A, rhs, M{j});
	  %      elapsed{j}				= elapsed{ j } + cputime - tstart{ j };

	  %  end

  end

  result(:,:,1)	= reshape(u{1}, dims);
  result(:,:,2)	= reshape(u{2}, dims);
  result(:,:,3)	= reshape(u{3}, dims);

	  
end
