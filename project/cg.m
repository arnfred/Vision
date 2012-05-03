% Assuming the P and R are all zero
% Translation of variables:
%
% X = B
% R = Diagonal
% P = L
% B = K
% b = rhs
%
% Original notation
% Computes c = A\b with A = X'*R*X + B'*P*B by Linear Conjugate Gradients (LCG) 

% My notation
% Solve system of: (B' * diagonal * B + P' * R * P) * u = rhs
%
% c is the solution, k is the amount of iterations, and res is the distance
function [c,k,res] = cg(X, R, b, epsilon)


	% The values used in the normal evaluation
	kmax	= 500;
	c		= zeros(size(b));
	tol		= 1e-3;
	prec	= [];


	% Setting up variables
	nb 		= norm(b,'fro');                        % |b|
	B		= matDiag(zeros(size(b)));
	res 	= zeros(kmax+1,1);


	% Check if fast Jacobi-DFT preconditioner is possible.
	% If it is, then we use as preconditioner the circulant matrix
	% M = mean(R)*X'*X + mean(P)*B'*B
	% Then M\r is implemented as pointwise division in the DFT domain.
	dX 		= diagFAtAFt(X,'cheap'); 
    R_mean 	= mean(R);
	d 		= R_mean*dX + epsilon;
	F 		= matFFTNmask(true(size(d))); d = d(:); % DFT matrix
	mvmMinvfun = @(r) mvmMinv(F,d, r);


	% Prepare fast function for matrix multiplication
	mvmAfun = @(c) mvmA(X, R, B, epsilon, c);
	r 		= b - mvmAfun(c);                       % r = b-A*c;
	z  		= r;
	res(1) 	= norm(r,'fro')^2;

	z = mvmMinvfun(r);
	rs = r(:)'*z(:);

	for k=1:kmax
		w = mvmAfun(z);
		alpha = rs/(z(:)'*w(:));
		c = c + alpha*z;
		r = r - alpha*w;
		res(k+1) = norm(r,'fro')^2;
		s = mvmMinvfun(r);
		rs_old = rs;
		rs = r(:)'*s(:);
		beta = rs/rs_old;
		z = s + beta*z;
		if sqrt(res(k+1))<tol*nb
			break
		end
	end
  
	res = sqrt(res(1:k+1));

	function v = mvmA(X,R,B,P, w)             % MVM with A where A = X'*R*X + B'*P*B
		Xw = X*w; 
		if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
		Bw = B*w; 
		if numel(P)==numel(Bw), PBw = P(:).*Bw; else PBw = P*Bw; end
		v = [X']*RXw + [B']*PBw;     % the [] around the transpose is needed by Octave
	end

	function z = mvmMinv(F,d, r)                             % M\r in the DFT domain
		if numel(r)>numel(d)
			z = cx2re([F']*((F*re2cx(r))./d));
		else
	  		z = [F']*((F*r)./d);
  		end
	end
end
