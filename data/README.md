# FDCVV Benchmark and Finite Element Data

This directory contains benchmark, ranking, computational-cost and finite
element datasets for FDCVV and comparison methods. Start with
`file_index.csv` for relative filenames, row counts, keys, figure mappings and
SHA256 hashes; `data_dictionary.csv` lists every column, storage type, unit and
observed missing count. `provenance.csv` binds every packaged table to its source.
Raw tables retain the supplied exports. Derived FE summaries and intervals are
recomputed as described below; provenance.csv binds each table to its source.

## Reading and structure

Files use UTF-8 (some source files include a UTF-8 BOM), comma delimiters, one
header row, decimal points, and standard CSV quoting. Scientific notation is
accepted. Empty fields or NaN denote missing values, never zero. Boolean fields
retain their original True/False or 0/1 representation. Keys are listed in the
index, separated by `|`. Data directories use relative paths and work after
unzipping on another computer. Absolute source paths retained in provenance
columns describe the original archive and are not import dependencies.

Python example: `pd.read_csv('numerical/baseline_endpoints540.csv')`.
MATLAB example: `readtable('finite_element/validation10.csv')`.

## Experimental units and identifiers

The numerical directories use GL1D, Peaks2D and H3_historical for Gramacy-Lee,
Peaks and the historical response H3(2*u-1), respectively. The ranking archive
also retains the original long-form benchmark identifiers. These are the same
three benchmark families. Each repeat denotes a complete trajectory. FDCVV
accuracy repeats are 106-135; other baseline methods use 1-30. Ablation and
lambda comparisons use matched trajectories. Steps and model sizes are not
additional independent repeats. The 75 LOOCV models comprise five sizes on
five trajectories for each benchmark. K_retention_figure.csv contains the
percentages plotted in Fig. 5; C1 retains the reference maximum-score region,
and C2 counts retained reference top-two positions. Panel a pools the 75 model
instances and panel b uses the 15 final-size models.

The frozen full configuration uses K=4,3,2 with transitions at enrichment
fractions 1/3 and 1/2, lambda0=0.35 and p_cross=0.10, and within-region
normalization. Version/config labels identify each control or lambda setting.

## Quantities and statistical suffixes

Numerical `common_*` fields use the fixed common test set. RMSE has function
response units. Numerical NRMSE is RMSE divided by the fitted model's predicted
response range on that test set. R2 is 1-SSE/SST. Native fields preserve original
exports but do not supply the adopted common-test comparisons.

The sensitivity gap is 100*(configuration mean / benchmark grid-minimum mean - 1).
The maximum and arithmetic mean of the three gaps are worst_gap_pct and
mean_gap_pct. Each metric has its own minimum. Descriptive rank ranks these
summaries. In paired ablation tables, difference is control minus Full and
relative_pct is 100*(control mean/Full mean-1). Wins mean lower Full NRMSE.
The Holm suffix adjusts across the 21 ablation comparisons. CI low/high fields
are pointwise 95% bootstrap bounds; see the Numerical supplement for resampling.
Mean, median, min and max suffixes describe the named quantity; std is sample
standard deviation. Summary n is the number of trajectories in that group.
LOOCV exact denotes the saved reconstruction reference; fast denotes the matrix
score. Model-size rows are repeated observations of the same trajectory.

Stopping decisions were stored before evaluation on the common test set.
Outcome fields compare stop and full-budget models; saved_points counts avoided
enrichments. triggered indicates rule satisfaction, including at the full budget.

## Complete computational cost

`timing/verified_trajectories.csv` contains 27 complete measurements, three per
method and benchmark. local repeat 1-3 maps to source_repeat 106-108 for FDCVV
and 1-3 for LC2015/CVVor. These are complete trajectory measurements, not 30
independent step repeats. Both initial_fit_s and refit_s include hyperparameter
optimization. core_s includes initial_fit_s + selection_s + truth_s + refit_s
and internal transition time. The latter is the residual after subtracting
the four measured components from core_s.
wait_s, test_s and trace_io_s are recorded separately. load_attempts.csv has
928 formal attempts; accepted=1 selects 837 adopted operations, while 91
attempts are excluded from the adopted cost. Nine warmups are excluded.

## Finite element data: corrected trajectory, 22 October 2026

The author supplied a new trajectory recalculated from training point 13 onward
and confirmed that these responses were obtained from manual finite element
calculations. Rows 1-12 are retained; the new complete training set has 30
positive responses. training_response_provenance30.csv records this origin.
Old processing labels are retained in the historical local audit, not attached
to the changed inputs of this replacement trajectory.

Training and validation are separate sets of 30 and 10 points. TD is delay in
seconds; OMEGA2 is propellant mass in kg; MH is piston mass in kg (M_C in the
manuscript). Responses are in metres, _um fields in micrometres, and prediction
variance in square metres. validation10.csv contains the final model predictions;
stage_validation60.csv contains six stages (5,10,15,20,25,30) on the same ten
validation inputs. stage_metrics26.csv retains all stages 5-30; the six stages
with pointwise predictions and final_metrics1.csv are recomputed from those
pairs. Other stages retain the supplied aggregate results. Missing MAE/max-error
values at the other twenty stages are blank, not zero.

FE NRMSE divides RMSE by the fixed observed validation-response range. The final
R2 is 0.9322546102, RMSE 11.90110528 um and NRMSE 8.503521106 percent.
bootstrap_intervals3.csv uses 100000 paired resamples of the ten validation pairs,
NumPy default_rng PCG64, seed 26100506, with the model and original response
range fixed. The 2.5th and 97.5th percentiles use linear interpolation.
bootstrap_protocol.json records the complete protocol and prediction-table hash.

The 30-point model was refitted with the frozen MS fitting procedure and seed
261016030. Its ten predictions agree with the supplied table within
1.24e-17 m. model_binding.json records the exact input and prediction hashes.

## Response-evaluation timing

Table 9 uses the same ten independent validation inputs. Their archived finite
element status-file intervals total 147810.361884 seconds (41.05843 h). Internal
solver wall-clock records total 147606 seconds and are retained separately.
The fitted new model was timed in single-threaded MATLAB R2022b: 20 warm-up
sweeps followed by 30 technical blocks of 100 ten-input sweeps. The median
per ten inputs is 0.730953 ms; fitting took 4.9934929 s separately. Timing
includes input normalization and mean/variance prediction, excluding model
fitting, loading, warm-up and output. Raw blocks and environment are included.

training_runtime30.csv preserves the supplied adaptive-duration export. Its
numeric durations match the preceding export and have not been independently
rebound to the recalculated training trajectory. This table is not used for
Table 9. The initial five duration fields are missing.

## Scope

The forty CSV data tables retain all numerical, ranking and complete-loop timing
evidence. Only the FE trajectory, derived summaries, prediction measurements and
related metadata were replaced. The manuscript and supplements describe the
scientific methods. Absolute source paths document local provenance and are
not required to read the portable tables.
