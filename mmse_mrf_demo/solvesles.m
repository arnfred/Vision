
% Matlab script for solving a bunch of linear systems based on the fields of
% experts mrf
function [result elapsed] = solvesles()

  % Dimensions
  dims = [82 82];

  % Add path
  addpath('mat');

  % Get mrf
  mrf = learned_models.cvpr_3x3_foe;

  % Express B in terms of the filters of the mrf
  B = ones(0,prod(dims));

  % For each filter, create a fake convolution matrix
  for i=1:mrf.nfilters
	  Bi = matConv2(mrf.filter(i), dims, 'valid');
	  B = vertcat(B,Bi); % B is 51200 x 6724
  end

  % Load Z matrices
  load('z.mat');

  % Now we should have a cell array called zes in scope
  [~, number_of_zs] = size(zes);

  % Set up timer
  elapsed = 0;
  tstart = cputime;

  % Set up return matrix as the last 4 matrices:
  result = zeros(dims(1), dims(2), 4);

  % Number of rows in diagonal (dims - 2 because the 'valid' fft transform gives back a smaller matrix) 
  N = prod(dims - 2) * mrf.nfilters; % size of filter (80*80) times number of filters (8)

  % Iterate over all Z's, solving the systems one by one
  for i = 1:number_of_zs
	  
	  % Get diagonal matrix
	  origDiagonal	= full(cell2mat(zes(i)));

	  % Cut off identity matrix saved at the end
	  diagonal		= matDiag(origDiagonal(1:N));
	  sqrtDiagonal	= matDiag(sqrt(origDiagonal(1:N)));
     
	  % Calculate right hand side
	  rhs			= B' * sqrtDiagonal * randn(N, 1); % n0 in notes

	  % Solve system of: (B' * diagonal * B + P' * R * P) * u = rhs
	  % Input are: B, diagonal, P, R, u, max MVM's, conditioner (I think), treshold, use_preconditioner
	  % These aren't exactly the same names as used in the documentation
	  P				= matDiag(zeros(1,prod(dims)));
	  R				= zeros(prod(dims),1);
	  maxMVMs		= 400;
	  conditioner	= zeros(prod(dims), 1);
	  treshold		= 1e-6;
	  use_precond	= 0;

	  % Plug in arguments and calculate away, time the operation
	  tstart		= cputime;
	  u				= linsolve_lcg(B, diagonal, P, R, rhs, maxMVMs, conditioner, treshold, use_precond);

	  % Calculate elapsed time
	  elapsed		= elapsed + (tstart - cputime);

	  % Save image vector in result matrix
	  n				= mod(i, 4) + 1;
	  result(:,:,n)	= reshape(u, dims);

  end
	  
end
