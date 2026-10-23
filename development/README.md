# Configuration development and trajectory reuse

The records in this directory document the development sequence for the
retained FDCVV configuration. `stage_index.csv` gives the configurations, run
identifiers, criteria, source records and later reuse. `source_index.csv` binds
each archived source record to its original and distributed hashes. Files in
`records` preserve their original contents, including historical labels and
source paths; those paths are not runtime dependencies.

## Earlier weight screening

The initial grid used lambda0 = 0.15, 0.20, 0.25, 0.30, 0.35 and p_cross =
0.10, 0.15, 0.20, 0.25, 0.30. Each setting used run identifiers 101-103 on each
of three benchmarks, giving 225 trajectories. K followed 4,3,2 with switches
at completed enrichment fractions 1/3 and 2/3; normalization was within each cell.

For theta=(lambda0,p_cross) and theta_ref=(0.25,0.20), the recorded criterion was

```text
mean_E[f,theta] = mean(final common-test RMSE over runs 101,102,103)
ratio[f,theta]  = mean_E[f,theta] / mean_E[f,theta_ref]
S[theta]       = (ratio[GL1D,theta] * ratio[Peaks2D,theta]
                  * ratio[H3_historical,theta])**(1/3)
```

The smallest S selected (0.35,0.10), with S=0.9053408729953407. Ties were
ordered by the largest benchmark ratio and then configuration name. This is
the ratio of within-benchmark means, followed by a cross-benchmark geometric
mean. `lambda_screening_endpoints225.csv`, `lambda_screening_ranking25.csv` and
`lambda_selection.json` contain the underlying results and recorded decision.

The selected weights, original weights and fixed lambda=0.5 were then compared
on identifiers 106-110. The 45 records are in `lambda_followup_endpoints45.csv`.
The historical JSON calls these identifiers `validation_repeats`; the same
identifiers subsequently participated in K-schedule development.

## Subsequent regional schedules

With weights (0.35,0.10) retained, six K432 alternatives were evaluated on
identifiers 106-110, alongside the original 1/3,2/3 schedule and archived controls.
Their switching fractions were (0.20,0.40), (0.20,2/3), (1/3,0.50), (1/3,5/6),
(0.50,2/3) and (0.50,0.75). The stage's `Full` label denotes the original
1/3,2/3 schedule; `K432_33_50` denotes the candidate that became the retained Full.

The 1/3,1/2 candidate was extended from five to fifteen runs (106-120).
Four neighbouring schedules, with switches after 8/15, 12/15, 10/13 and 10/17
completed additions, were also compared using identifiers 106-110. These are
counts out of a 30-addition budget, as recorded in `K_schedule_neighbours.json`.
The retained candidate and three weight controls were then extended to 30 runs
(106-135). Their 360 endpoints are included in the development records.

The K schedule was retained through these developmental comparisons of
benchmark errors and paired differences. The archive does not record a single
prespecified scalar selection rule for K. The numerical screening score S above
applies to the earlier weight selection, not to the subsequent K decision.
`configuration_freeze.json` records the retained 4,3,2 schedule with switches at
1/3 and 1/2 and weights (0.35,0.10).

## Later sensitivity and reuse

The later 25-setting grid uses the adopted 1/3,1/2 K schedule, lambda0 =
0.05,0.15,0.25,0.35,0.45 and p_cross = 1/30,0.10,0.20,1/3,0.50.
There are 30 trajectories per setting and benchmark. The maximum relative gap G
uses each benchmark's minimum mean RMSE in this later grid. It is a retrospective
sensitivity summary; it was not the criterion for the earlier weight screening.

The same 90 retained trajectories are used as FDCVV in the baseline, Full in the
ablation and the retained parameter cell in the grid. `trajectory_reuse90.csv`
lists the benchmark, run identifier, three dataset labels and endpoint values.
The grid marks these rows `reused=True`. These uses do not add independent runs.
All three development benchmarks and their common evaluation sets were reused
across these stages. No separate confirmation experiment was performed after
the complete configuration was fixed.

Archived `native_nrmse`, `nrmse` and `common_nrmse` fields retain their historical
definitions. The developmental RMSE records and screening ratios remain unchanged.
For the current numerical comparisons, use `true_nrmse` in `data/numerical`.
`python CHECK_DEVELOPMENT.py` checks the recorded weight ranking and trajectory
reuse without fitting models.
