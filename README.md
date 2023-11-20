# Issa2023
Code related to Issa et al - Nat Neuro 2023

### findbn.m
Given fluoresence trace y and neuropil signal neu, attempts to recover signals by assuming y = spk_est * h + B * coeff, where h is a kernel defined by time constants tau, B is the basis functions, and coeff is the associated constants. Thus spk_est*h is the estimate of F(t) and B * coeff is the estimate of the baseline. Algorithm iterates between solving for coeff (using lsqlin) and updating spk_est (using oasisAR2).

For convenience, we have included oasisAR2.m, oasisAR1.m, and thresholded_oasisAR2.m from https://github.com/flatironinstitute/CaImAn-MATLAB/tree/master/deconvolution/oasis. These files are needed for findbn to run. However, different spike inference methods can be used in place of Oasis if desired by replacing the appropriate code in the est_spk subroutine in findbn.m.

Please cite/reference Issa et al - Nat Neuro 2023.
If Oasis is used, please cite/reference Friedrich et al - PLOS Comp Bio 2017 and/or Giovannucci et al - eLife 2019.
