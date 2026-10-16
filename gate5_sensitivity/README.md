# Sensitivity and uncertainty

Supports Model page section 5.1 and supplementary S12. One self-contained
script carrying its own re-implementation of Models 3 to 7, so the analysis
does not depend on the Gate 4 scripts being importable.

```bash
pip install numpy scipy SALib
python3 nishvasa_sensitivity_simulation.py verify    # baseline check, seconds
python3 nishvasa_sensitivity_simulation.py partA     # tornado, Morris, Sobol
python3 nishvasa_sensitivity_simulation.py partB     # Monte Carlo
python3 nishvasa_sensitivity_simulation.py all
```

## Run verify first

A sensitivity analysis of a model that does not reproduce its own baseline
measures nothing. `verify` checks the independent re-implementation against
the modelling document and cross-checks the channel ODE solver:

| Quantity | Document | This code |
| --- | --- | --- |
| Removal at 20 ppm | 585 | 585.0 umol/min |
| Removal at 50 ppm | 715 | 712.7 umol/min |
| Removal at 100 ppm | 765 | 762.7 umol/min |
| Batch, 500 birds, 1 unit | 6.0 h | 6.04 h |
| Batch, 1500 birds, 3 units | 7.3 h | 7.34 h |
| Batch, 3000 birds, 6 units | 9.6 h | 9.59 h |

The 2 umol/min difference at 50 ppm is the coating effectiveness factor: this
code computes 0.9848 from the diffusion equation where the document rounds to
0.98. RK4 against LSODA on the channel integration agrees to 1e-10.

## What each part does

**Part A** -- local elasticities at the baseline, Morris screening over 50
trajectories, and Sobol decomposition at 2,048 base samples with 1,000
bootstrap resamples, plus a convergence replicate at 1,024. Sixteen inputs,
two outputs: removal rate at 50 ppm and batch time for the 1500-bird farm.

Batch time is decomposed on a log scale. The model is a product of factors,
so a log transform turns it into a sum and the variance shares total 100 %.

**Part B** -- 10,000 Latin hypercube draws, converted to a success
probability per unit count and the count needed for 90 % confidence. Writes
`S2_deployment_table.csv`.

## Reading the result

Four inputs -- enzyme loading, retained activity, specific activity and
coating thickness -- multiply together in the V_max,eff expression and
account for 68 % of the batch-time variance. They are one quantity: active
enzyme per square metre of wall.

Sobol shares scale with the ranges tested, so loading ranks above retention
partly because its range is fourfold against retention's 2.3-fold. The
elasticities, which do not depend on ranges, place the four within 3 % of
each other. Use both.

## Scope

The device layer only, Models 6 and 7. The achieved enzyme ratio, phosphate
inhibition of PPK2 and operating pH sit in Gate 3 and are not propagated
into it; S12.7 states why each is absent.
