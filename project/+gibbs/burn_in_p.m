function p = burn_in_p(mrf, u, treshold)

	% Get the energy of the last half of the samples
	u_energy	= mrf.energy(u(ceil(end/2):end));

	% Compute the epsr of the last half of the samples
	epsr		= pml.support.epsr(u_energy);

	% Check if the epsr is above the treshold
	p 			= (epsr > treshold);
end
