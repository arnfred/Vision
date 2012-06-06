function psnr = get_psnr(original, noisy, border, dims)

	ridx = 1+border:dims(1)-border; 
	cidx = 1+border:dims(2)-border;
	[rs, cs] = ndgrid(ridx, cidx);

	% indices of interior pixels
	ind_int = sub2ind(dims, rs(:), cs(:));

	% Get PS
    psnr = pml.image_proc.psnr(reshape(noisy(ind_int), size(original)), original);

end
