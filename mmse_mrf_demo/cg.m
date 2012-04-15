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
function [c,k,res] = cg(B, Diagonal, rhs, dims)

% The values used in the normal evaluation
kmax	= 400;
c		= zeros(size(rhs));
tol		= 1e-6;
prec	= [];

% Setting up variables
nb 		= norm(rhs,'fro');                        % |b|
K		= matDiag(zeros(1,prod(dims)));
L		= zeros(prod(dims),1);
res 	= zeros(kmax+1,1);

mvmAfun = @(c) mvmA(B, Diagonal, K, L, c);
r 		= rhs - mvmAfun(c);                       % r = b-A*c;
z  		= r;
res(1) 	= norm(r,'fro')^2;

for k=1:kmax
		w = mvmAfun(z);
		alpha = res(k)/(z(:)'*w(:));
		c = c + alpha*z;
		r = r - alpha*w;
		res(k+1) = norm(r,'fro')^2;
		beta = res(k+1)/res(k);
		z = r + beta*z;
		if sqrt(res(k+1))<tol*nb
			break
	end
end

function v = mvmA(X,R,B,P, w)             % MVM with A where A = X'*R*X + B'*P*B
  Xw = X*w; if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
  Bw = B*w; if numel(P)==numel(Bw), PBw = P(:).*Bw; else PBw = P*Bw; end
  v = [X']*RXw + [B']*PBw;     % the [] around the transpose is needed by Octave
