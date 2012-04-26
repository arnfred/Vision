
% Matlab script for solving a bunch of linear systems based on the fields of
% experts mrf
function [result elapsed] = solvesles()

  % Dimensions
  dims = [82 82];

  % Add path
  addpath('mat');
  addpath('inf');

  % Get mrf
  mrf = learned_models.cvpr_3x3_foe;

  % Express B in terms of the filters of the mrf
  % B = ones(0,prod(dims));
  B = mat([0 prod(dims)], 'valid');

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
  elapsed = zeros(3,1);
  tstart = zeros(3,1);

  % Set up return matrix as the last 4 matrices:
  result = zeros(dims(1), dims(2), 3);

  % Number of rows in diagonal (dims - 2 because the 'valid' fft transform gives back a smaller matrix) 
  N = prod(dims - 2) * mrf.nfilters; % size of filter (80*80) times number of filters (8)

  % Iterate over all Z's, solving the systems one by one
  % for i = 1:number_of_zs
  for i = 1:5
	  
	  % Get diagonal matrix
	  origDiagonal	= full(cell2mat(zes(i)));

	  % Cut off identity matrix saved at the end
	  diagonal		= matDiag(origDiagonal(1:N));
	  sqrtDiagonal	= matDiag(sqrt(origDiagonal(1:N)));
     
	  % Calculate right hand side
	  rhs			= B' * sqrtDiagonal * randn(N, 1); % n0 in notes

	  % Prepare preconditioners
	  A				= B'*diagonal*B;
	  diag_A		= diag(A);
	  d 			= mean(diag(diagonal)) * diagFAtAFt(B,'cheap');
	  F 			= matFFTNmask(true(size(d))); % DFT matrix
	  M1			= @(r) r;
	  M2			= @(r) r ./ diag_A;
	  M3			= @(r) mvmMinv(F, d(:), r);


	  % Plug in arguments and calculate away, time the operation
	  tstart(1)		= cputime;
	  [u1, iter1]	= conjugate_gradients(A, rhs, M1);
	  elapsed(1)	= elapsed(1) + (cputime - tstart(1));

	  % Now do the same, but with diag(A) preconditioning
	  tstart(2)		= cputime;
	  [u2, iter2]	= conjugate_gradients(A, rhs, M2);
	  elapsed(2)	= elapsed(2) + (cputime - tstart(2));

	  % Now with diagonal in fuirier space as preconditioning
	  tstart(3)		= cputime;
	  [u3, iter3]	= conjugate_gradients(A, rhs, M3);
	  elapsed(3)	= elapsed(3) + (cputime - tstart(3));

  end

  result(:,:,1)	= reshape(u1, dims);
  result(:,:,2)	= reshape(u2, dims);
  result(:,:,3)	= reshape(u3, dims);

	  
end

  % M\r in the DFT domain (From hannes nickish code)
  function z = mvmMinv(F,d, r)
	  if numel(r)>numel(d)
		  z = cx2re([F']*((F*re2cx(r))./d));
	  else
		  z = [F']*((F*r)./d);
	  end
  end
