function u = test_sampler()

  mrf = learned_models.cvpr_3x3_foe;
  dims = [82 82];
  mrf.imdims = dims% + [2 2];
  mrf.update_filter_matrices = true;
  mrf.conv_method = 'circular';

  u = gibbs.sample(mrf, dims);
end
