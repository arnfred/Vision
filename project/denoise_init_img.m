function [img_noisy img_clean] = denoise_init_img(sigma)

	% Set some variables
	scale	= 1/8; % As used in Darmstadt paper
	img_names = {'barbara.png', 'boat.png', 'fingerprint.png', 'house.png', 'lena.png', 'peppers256.png'};
	idx = 5; % choose index of img_names

	% Load image and scale
	img_clean = double(imread(sprintf('images/denoising/%s', img_names{idx})));
	img_clean = imresize(img_clean, scale);

	% Create noisy test_image
	img_noisy = img_clean + sigma * randn(size(img_clean));
end
