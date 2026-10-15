# Sensitivity and uncertainty

To be added. The analysis is reported in Model page section 5.1 and
supplementary S12. Four methods on two outputs:

| Method | Scope |
| --- | --- |
| Tornado | +/- 20 % on each of sixteen inputs, local elasticities |
| Morris | 50 trajectories, screening |
| Sobol | 2,048 base samples, 1,000 bootstrap resamples |
| Monte Carlo | 10,000 Latin hypercube draws |

Outputs: removal rate at 50 ppm for one unit, and batch time for the
1500-bird farm. Batch time is analysed on a log scale, because the model
is a product of factors and the variance shares then sum to 100 %.

The code should reproduce the baseline check in S12.1 before anything
else: removal 713 against 715 umol/min, all three batch times exact.
