% Removes the border around an image and reshapes it
function u = get_img(u_padded, border, dims)
	ridx = 1+border:dims(1)-border; 
	cidx = 1+border:dims(2)-border;
	[rs, cs] = ndgrid(ridx, cidx);

	% indices of interior pixels
	ind_int = sub2ind(dims, rs(:), cs(:));
	dims
	border

    u = reshape(u_padded(ind_int), dims - 2*border)/255;
end
