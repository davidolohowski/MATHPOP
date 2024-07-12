# Data

Raw data and post-inference data from the original MATHPOP paper.

* `fake_data`: Artificial star test data used to obtain photometric measurement uncertainties and completeness fraction of images under DOLPHOT for the PIPER survey.
* `GC_prob`: Combined results from fitting non-parametric finite-mixture model to the DOLPHOT point source data to obtain the probabilistic GC catalog.
* `point_source_data`: Point source data within the PIPER survey obtained from SExtractor (`PS_dat_Jans.csv`) and DOLPHOT (all other files).
* `prob_GC_data`: Probabilistic GC catalogs: files containing `_Jans` are catalogs based on SExtractor, others are based on DOLPHOT.
* `sim`: Summarized data results from simulation study.
* `summary_results.RDS`: Summarized data results fitted to the real data.
* `v###`: Every file contains the MCMC results fitted to the SExtractor (with `_Jans`) and DOLPHOT (without `_Jans`) GC data.