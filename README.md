## Model overview

![Model architecture](results/figures/model_architecture_polished.png)

The workflow predicts transcriptional responses to perturbations that
were excluded from training. Perturbation identities are split before
pseudobulk construction to prevent target leakage.

## Project summary

![Portfolio summary](results/figures/portfolio_summary_polished.png)

## Biology-informed model comparison

![Biology-informed comparison](results/figures/biology_prior_model_comparison_polished.png)

The lowest mean MAE in this experiment was obtained by
**Baseline**
(**0.0312**).

## Dose-aware model comparison

![Dose-aware comparison](results/figures/dose_model_comparison_polished.png)

The lowest mean MAE in this experiment was obtained by
**Hill gradient**
(**0.0358**).

## Predicted-versus-observed response magnitude

![Magnitude calibration](results/figures/predicted_vs_observed_magnitude_by_model.png)

Points close to the diagonal indicate that the model predicts the
overall transcriptomic response magnitude accurately. Points below the
diagonal indicate systematic underprediction.

![Calibration error](results/figures/magnitude_calibration_error_by_model.png)

## Dose-response trajectories

![Dose trajectories](results/figures/dose_response_trajectories_all_models.png)

## Generalization across held-out targets

![Held-out target performance](results/figures/heldout_target_performance_polished.png)

> Biology-informed and dose-aware experiments are shown separately
> because raw MAE values should only be directly compared when the
> underlying test data and pseudobulk construction are identical.
