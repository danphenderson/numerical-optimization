### Overview
--------------

Interested in the critical structure of smooth functions $f: \R^n \rightarrow \R$, e.g. consider
We seek a method of determining which functions are rich in saddle points. We belive that having
large negative curvature near a saddle incereases likelyhood of tacking for Quasi-Newton methods
that approximate the region with a quadratic Taylor Polynomial; we attempt to characterize our set
of test functions to better understand their behavior

### Numerical Experiment

**Implementation Details**
We build symbolic expressions for $\nabla f$ and $\nabla^2f$ in Mathematica:
```Mathematica
BuildGradAndHess[f_, n_] := Module[
   {x, i, var, df, ddf},
   var = Array[x, n];
   df = D[f[var], {var}];
   ddf = D[f[var], {var, 2}];
   Quiet[
    {Function[Evaluate[x], Evaluate[df /. x[i_] :> Part[x, i]]],
     Function[Evaluate[x], Evaluate[ddf /. x[i_] :> Part[x, i]]]}
    ]
   ];
```

We generate a sample box around our initial iterate, $x_0$, to define our sampling space:
$$ \mathbb{B} = \{ \vec{x} \in \mathbb{R}^n : (\vec{x_o})_i - R_i \leq x_i \leq (\vec{x_o})_i + R_i  \}$$
where $R_i = 2 |(\vec{x_o})_i|$ provided that $|(\vec{x_o})_i| > tol$, or, $R_i = 1$. Observe
```Mathematica
SampleFromBox[x0_, s_Integer] := 
 Module[{n = Length[x0], Tol = 0.1, Scale},
  Scale = 2.0*Abs[x0];
  Do[If[Scale[[i]] <= 2.0*Tol, Scale[[i]] = 1], {i, 1, n}];
  Table[
   RandomReal[{-1, 1}, n]*Scale + x0,
   {s}]
]
```

We declare the utility functions
```Mathematica
CalculateIndex[eigenvalues_] := Module[
  {signs = Sign[eigenvalues]},
  {
    Count[signs, 1],  (* Positive eigenvalues *)
    Count[signs, 0],  (* Zero eigenvalues *)
    Count[signs, -1]  (* Negative eigenvalues *)
  }
]
```

Now we modularize our main routine as:

```Mathematica
AnalyzeFunctionCriticalPoints[f_, z_, s_Integer, opts : OptionsPattern[]] := Module[
  {n, gradFunc, hessFunc, samples, sampledEvals, indices, tallyIndices, plotOption},
  
  n = Length[z];
  {gradFunc, hessFunc} = BuildGradAndHess[f, n];
  samples = SampleFromBox[z, s, OptionValue["Tolerance"], OptionValue["ScaleFactor"]];
  
  sampledEvals = Eigenvalues[hessFunc /@ samples];
  indices = CalculateIndex /@ sampledEvals;
  tallyIndices = Sort[Tally[indices]];
  
  (* Plotting if enabled *)
  If[OptionValue["PlotEigenvalues"],
    ListPlot[sampledEvals, PlotLabel -> "Eigenvalues at Sampled Points"],
    Null
  ];
  
  (* Return results *)
  <|
    "TallyIndices" -> tallyIndices,
    "SampledEvals" -> sampledEvals
  |>
]
```






### Things to Consider
- TODO: use the idea of a determinant to measure the 'volume' of our open ball we are sampling?
  - How would this relate to $\text{det} \nabla^2 f(\vec{s_i})$ at each sample point.
  - Also, consider how this relates to the Trust Region because we are ulitmately trying to understand to taking nature that occurs in the next few steps. Well, our step criterion may assert that we can only move so far away in any finite amount of steps. So really, we may be insterested in learning more about 'tacking' after we find some good canidates.
- TODO: Introduce mathematical notation:
  - define sample space, measure of sample space (or volume from det) 
- TODO: We ultimately must classify $f$'s structure statistically..
  - Create a criterion, include various metrics to assess topology of
    sample space.
  - Consider how we may classify $\nabla^2 f(\vec{s_i})$ sadles which have dominant negative-curvature directions.
    - Note, we may never actually get to saddles.
- TODO: Try to see if insight is gained from looking at one problem in variying dimension.