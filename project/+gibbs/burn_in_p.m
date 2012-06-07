function p = burn_in_p(mrf, energy, treshold)

	% The number of iterations
	iter		= numel(energy);

	% If we have iterated less than 10 times, we continue
	if (iter < 10) p = 1; return; end
		
	% Compute the epsr of the last half of the samples
	% epsr		= pml.support.epsr(energy(end-9:end))
	variance	= var(energy(iter-9:iter))

	% Check if the epsr is above the treshold
	% p 			= (epsr > treshold);
	p			= (variance > treshold);
end
