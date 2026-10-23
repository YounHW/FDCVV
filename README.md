# FDCVV: Adaptive Sampling with Fast LOOCV and Voronoi Partitioning

This repository provides the FDCVV adaptive sampling implementation, numerical
comparison methods, fixed test sets, and machine-readable benchmark and
finite element results.
Download the pinned upstream dependency first using the instructions below.
Use this folder as the MATLAB current directory.

## Algorithm configuration

FDCVV retains Fast LOOCV and Voronoi screening. For enrichment step `m=1:30`,
`q=(m-1)/30` determines K: 4 for steps 1–10, 3 for 11–15, and 2 for 16–30.
The variance weight is `lambda=.35+.65*(m/30)^alpha`, where
`alpha=log(.15/.65)/log(.10)`. Distance has weight `1-lambda`.
Both terms are normalized within each selected Voronoi cell. Hyperparameters
are re-estimated by the original multistart fitting procedure after each addition.

`matlab/final_fdcvv/FDCVV_function.m` is the current implementation.
`frozen_configuration.json` records the source configuration. The execution
statements preserve the frozen implementation; obsolete introductory comments
and working-directory-dependent `addpath` directives have been removed.
`source_bindings.csv` records original and packaged file hashes and these changes.

## Requirements

Validated with MATLAB R2022b on Windows. Full numerical runs use Optimization
Toolbox (`fmincon`), Global Optimization Toolbox (`MultiStart`) and Statistics
and Machine Learning Toolbox (`lhsdesign`, distance functions). A serial run
does not require Parallel Computing Toolbox. Version and toolbox details are
recorded in the validation reports. Sampling code is MATLAB. Python 3.11 with
NumPy, pandas and SciPy reconstructs the archived numerical statistics;
`requirements-statistics.txt` gives the tested versions. The checksum verifier
uses only Python's standard library.

## Original framework dependency

The original framework is [FuhgJan/StateOfTheArtAdaptiveSampling](https://github.com/FuhgJan/StateOfTheArtAdaptiveSampling).
This release uses commit `f3d9b779ba3ac97f2cdd85f6cc93ef2a1e6a3626`. This is the dependency revision verified
for this package, not a claim about the historical ZIP download. Its unchanged
files are not redistributed here. `upstream_comparison.csv` records all 91
previous base files: 67 supplied by upstream and 24 local additions or changes.
Comparisons allow line endings, trailing whitespace and redundant relative
addpath directives. All numerical changes remain in the local overlay.

Download and extract this exact upstream snapshot:
https://github.com/FuhgJan/StateOfTheArtAdaptiveSampling/archive/f3d9b779ba3ac97f2cdd85f6cc93ef2a1e6a3626.zip

Or use Git:
```sh
git clone https://github.com/FuhgJan/StateOfTheArtAdaptiveSampling.git
git -C StateOfTheArtAdaptiveSampling checkout f3d9b779ba3ac97f2cdd85f6cc93ef2a1e6a3626
```

At the start of each MATLAB session, set the extracted upstream directory:
```matlab
setenv('FDCVV_UPSTREAM_ROOT','C:/your/path/StateOfTheArtAdaptiveSampling');
```
It must be the folder containing `src`, not `src` itself. SETUP_PACKAGE checks
the dependency files against `upstream_lock.json`, loads them, then prioritizes
the local modified files. Keep the downloaded upstream directory unmodified.
The original framework's relative addpath calls may warn when running from an
output directory; required folders are already loaded by absolute path.

The original framework's MIT notice is retained in
`THIRD_PARTY_LICENSE_Fuhg_MIT.txt` for modified third-party files. Cite Fuhg,
Fau and Nackenhorst, "State-of-the-Art and Comparative Review of Adaptive
Sampling Methods for Kriging" when using the framework. The local LC2015
directory supplies a single-fidelity implementation and its framework adapter.

## Quick verification

Open a fresh MATLAB session for each command (the two methods share a class
name and use path-selected implementations):

```matlab
CHECK_FDCVV_REPLAY
```

This checks 15 saved model states: three benchmarks at steps 10, 11, 15, 16 and
30 of source repeat 106. It checks selected coordinates, top-K regions, scores,
regional winners, the K schedule, lambda values, and restoration of the random
state. The fixture contains saved X, Y, theta, R and expected selection traces.
No hyperparameter optimization or new true-response evaluation is performed.

In a separate fresh MATLAB session:

```matlab
CHECK_FIXED_MODEL('LC2015')
```

This checks fixed-parameter analytical LOO against explicitly reconstructed
fixed-parameter models, candidate-score maximization, search convergence and
resolution of dependencies in this package and the pinned upstream directory. Reports appear under
`validation`. These checks test the portable implementation; the archived
30-run results remain in `data`.

## Run a complete numerical trajectory

In a fresh MATLAB session, run one of:

```matlab
RUN_ONE_CASE('FDCVV','GL1D',106)
RUN_ONE_CASE('LC2015','Peaks2D',1)
```

These are full runs with 30 added observations and hyperparameter optimization.
Supported methods are FDCVV, LC2015, CVVor, SSA, MaxVariance and LHS_rebuild.
Benchmarks are GL1D, Peaks2D and H3_historical. Initial sizes are 5, 10 and 15;
final sizes are 35, 40 and 45. Hartmann uses H3(2*u-1) on normalized input u.

Use FDCVV source IDs 106–135, and IDs 1–30 for the other methods.
`configurations.csv` provides each initial-design and fitting seed.
Initial designs are fixed within each benchmark. For FDCVV, benchmark index
b=1,2,3 and source ID r give fitting seed `281010000+b*100000+r` and
candidate seed `271010000+b*100000+r*100+m`. Candidate selection restores the
incoming RNG state. K progress uses `(m-1)/30`; lambda progress uses `m/30`.

Results are written to `results/<method>/<benchmark>/rep_<ID>`.
An existing result directory is never overwritten. Files include launch inputs,
model history, native convergence metrics, common-test convergence metrics and
an acceptance record. FDCVV also saves `trace_XX.mat` with the selected point,
regional winners, fast scores, K, lambda and RNG state. Common-test NRMSE uses
the fixed true-response range and is exported as `true_nrmse`.
The same denominator is used for every method and stage within a benchmark.
`common_nrmse` retains the legacy prediction-range definition, and
`native_convergence.csv` retains the original framework exports.
The FE tables use the observed validation-response range.
Saved data reports remain unchanged when a new trajectory is run.

## Reconstruct numerical statistics

Run the reporting-metric check in MATLAB; this does not fit a model:

```matlab
CHECK_COMMON_METRICS
```

Reconstruct the numerical tables from the saved run-level RMSE and fixed
test responses:

```sh
python -m pip install -r requirements-statistics.txt
python REBUILD_STATISTICS.py --check
```

Outputs are written to `results/rebuilt_statistics`, leaving the supplied data
unchanged. They include baseline endpoint and convergence summaries, ablation
effects and Holm corrections, paired-bootstrap intervals, the weight grid and
stopping summaries. The check compares all nine affected numerical tables with
the supplied exports. No model fitting or new function evaluations are needed.
Resampling seeds, draw counts, input hashes and software versions are recorded
in the generated `verification.json`. The statistical unit is a complete
trajectory, conditional on the fixed initial design. Technical timing repeats
are reported separately.

The primary numerical field is `true_nrmse = common_rmse / (max(truth)-min(truth))`.
It is a dimensionless fraction; multiply by 100 to report percent.
`data/normalization_protocol.json` records the fixed denominators and test-set
hashes. The legacy `common_nrmse` and `native_nrmse` fields remain available for
reading historical exports and are not used as the primary numerical metric.

## Configuration development and shared trajectories

`development/README.md` explains the earlier weight screening, subsequent K
schedule comparisons and later sensitivity grid. `development/stage_index.csv`
links each stage to its configurations, run identifiers, criterion and archived
records. The earlier weight criterion compares benchmark-wise mean RMSE with
the original `(0.25,0.20)` setting and takes their geometric mean.

The same 90 FDCVV trajectories (30 per benchmark, identifiers 106-135) are used
in the six-method baseline, the Full ablation configuration and the retained
weight-grid cell. `development/trajectory_reuse90.csv` maps these shared records.
Some trajectories participated in parameter and schedule development. No separate
confirmation experiment was performed after the final configuration was fixed.
The later maximum-gap statistic describes sensitivity under the adopted schedule.

## Data and provenance

`data/README.md`, `data/file_index.csv` and `data/data_dictionary.csv` describe
40 CSV data tables, their units, keys and source mappings. `tests` contains
the fixed common numerical test sets and the compact FDCVV replay fixture.
The code package runs with relative paths; original absolute source paths in
provenance records identify the archive and are not runtime dependencies.
The FE component consists of the reported input/response, validation, interval
and runtime data; this numerical entry point does not launch Abaqus.

`SHA256SUMS.csv` lists distributed files. Run `python VERIFY_FILES.py` to check
them using Python's standard library. The MATLAB verification reports describe
checks actually executed for this package. Computational performance is
documented by the controlled measurements in `data/timing`.

## License and citation

FDCVV additions and modifications are distributed under the MIT License in
`LICENSE`. Modified upstream framework portions retain Jan N. Fuhg's copyright
and the original notice in `THIRD_PARTY_LICENSE_Fuhg_MIT.txt`.
`CITATION.cff` provides the software citation for this repository. Please also
cite the upstream framework and the methods used, as described above.
A release version and DOI can be added to the citation file when assigned.

## Finite element data update, 22 October 2026

The finite element dataset includes a 30-point training trajectory recalculated
from point 13 onward and ten independent validation inputs. The final model
has R2 = 0.93225, RMSE = 11.9011 um and NRMSE = 8.5035%.
The update includes all 26 stage summaries, predictions at six stages,
paired-bootstrap intervals, model/data binding records and prediction timings.
The median time to evaluate all ten validation inputs is 0.730953 ms;
model fitting took 4.9934929 s separately. The algorithm implementation retains
its verified version. Numerical reporting now uses the fixed true-response
normalization described above; raw RMSE, R2 and saved trajectories are unchanged.

See `data/README.md` for units, statistical definitions and timing scope.
Run `python VERIFY_FILES.py` after downloading or cloning to verify all
distributed files. `.gitattributes` preserves the exported bytes across platforms.
