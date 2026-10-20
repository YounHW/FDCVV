# CM-26-0575 Code and Data Supplement

This package supplies the adopted FDCVV implementation, numerical comparison
methods, fixed common test sets, and the machine-readable manuscript data.
Download the pinned upstream dependency first using the instructions below.
Use this folder as the MATLAB current directory.

## Adopted algorithm

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
recorded in the validation reports. Source code is MATLAB; Python is optional
for the standard-library checksum verifier.

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
directory supplies the single-fidelity method and adapter used in this study.

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
the predicted response range; the FE tables use the observed validation range.
Saved data reports remain unchanged when a new trajectory is run.

## Data and provenance

`data/README.md`, `data/file_index.csv` and `data/data_dictionary.csv` describe
36 CSV data tables, their units, keys and manuscript mappings. `tests` contains
the fixed common numerical test sets and the compact FDCVV replay fixture.
The code package runs with relative paths; original absolute source paths in
provenance records identify the archive and are not runtime dependencies.
The FE component consists of the reported input/response, validation, interval
and runtime data; this numerical entry point does not launch Abaqus.

`SHA256SUMS.csv` lists distributed files. Run `python VERIFY_FILES.py` to check
them using Python's standard library. The MATLAB verification reports describe
checks actually executed for this package. Computational performance is
documented by the separately controlled measurements in the timing supplement.
