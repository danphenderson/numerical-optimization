Implement `programs-probe.ipynb` in Julia, building on `mathematica/programs-probe.nb`
  - Target: Final list of test problems with dimensions, number of positive and negative eigenvalues, and median.
Implement `programs-solve.ipynb` in Julia, building on `mathematica/programs-solve.nb`
  - Target: Statistics on stalled minimizations, evaluations, and gradients.
Implement `scripts/adnlpmodels-to-mathematica.py` converting test set to Mathematica for symbolic analysis.
  - Target: Replace the current Mathematica test-set.



**Generate Report**

It should include the following:

*Table 1*
- random points, number of positive and negative eigenvalues and median.
  
*Table 2*
- stuff from stalled minimizations, evals(+/-), gradient.
- Quasi Newton scheme defaults in mathematica. 

Also, mention any conclusions from varying dimensions. Additionally,
randomly pick a few interesting programs to look into in depth.
E.g., reproduce the cyclic behvaior of the Generalized Rosenbrock function.