from __future__ import annotations

from pathlib import Path
import math
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

try:
    from IPython.display import display
except ImportError:
    def display(obj):
        print(obj)


# ============================================================
# 1. Resolve repository paths safely
# ============================================================
try:
    SCRIPT_PATH = Path(__file__).resolve()
    PROJECT_ROOT = SCRIPT_PATH.parents[1]
except NameError:
    CURRENT_DIR = Path.cwd().resolve()

    if CURRENT_DIR.name in {
        "notebooks",
        "scripts",
    }:
        PROJECT_ROOT = CURRENT_DIR.parent
    else:
        PROJECT_ROOT = CURRENT_DIR


SOURCE_OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "source_outputs"
)

RESULT_DIR = (
    PROJECT_ROOT
    / "results"
)

FIGURE_DIR = (
    RESULT_DIR
    / "figures"
)

FIGURE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


print("Project root:", PROJECT_ROOT)
print("Source outputs:", SOURCE_OUTPUT_DIR)
print("Figure directory:", FIGURE_DIR)


# ============================================================
# 2. Input paths
# ============================================================
INPUT_FILES = {
    "biology_summary": (
        SOURCE_OUTPUT_DIR
        / "multiheldout_model_summary.csv"
    ),
    "biology_per_target": (
        SOURCE_OUTPUT_DIR
        / "multiheldout_per_target_runs.csv"
    ),
    "biology_overall": (
        SOURCE_OUTPUT_DIR
        / "multiheldout_overall_runs.csv"
    ),
    "biology_paired": (
        SOURCE_OUTPUT_DIR
        / "multiheldout_paired_comparisons.csv"
    ),
    "biology_wins": (
        SOURCE_OUTPUT_DIR
        / "multiheldout_win_summary.csv"
    ),
    "gradient_summary": (
        SOURCE_OUTPUT_DIR
        / "gradient_model_summary.csv"
    ),
    "gradient_per_target": (
        SOURCE_OUTPUT_DIR
        / "gradient_per_target_runs.csv"
    ),
    "gradient_per_bin": (
        SOURCE_OUTPUT_DIR
        / "gradient_per_bin_runs.csv"
    ),
    "gradient_overall": (
        SOURCE_OUTPUT_DIR
        / "gradient_overall_runs.csv"
    ),
    "gradient_paired": (
        SOURCE_OUTPUT_DIR
        / "gradient_paired_comparisons.csv"
    ),
    "gradient_wins": (
        SOURCE_OUTPUT_DIR
        / "gradient_win_summary.csv"
    ),
}


REQUIRED_KEYS = [
    "biology_summary",
    "biology_per_target",
    "gradient_summary",
    "gradient_per_target",
    "gradient_per_bin",
]


print("\nInput-file check:")

for key, path in INPUT_FILES.items():
    print(
        f"{key}:",
        "FOUND"
        if path.exists()
        else "NOT FOUND",
        "—",
        path,
    )


missing_required = [
    key
    for key in REQUIRED_KEYS
    if not INPUT_FILES[
        key
    ].exists()
]

if missing_required:
    missing_text = "\n".join(
        f"- {key}: {INPUT_FILES[key]}"
        for key in missing_required
    )

    raise FileNotFoundError(
        "Required files are missing:\n"
        f"{missing_text}"
    )


# ============================================================
# 3. Loading and column helpers
# ============================================================
def read_result_table(
    path: Path,
) -> pd.DataFrame:
    table = pd.read_csv(
        path
    )

    unnamed = [
        column
        for column in table.columns
        if str(
            column
        ).startswith(
            "Unnamed:"
        )
    ]

    if (
        "model" not in table.columns
        and unnamed
    ):
        table = table.rename(
            columns={
                unnamed[0]: "model",
            }
        )

    return table


def first_existing_column(
    table: pd.DataFrame,
    candidates: list[str],
    table_name: str,
) -> str:
    for column in candidates:
        if column in table.columns:
            return column

    raise KeyError(
        f"{table_name} does not contain any of {candidates}.\n"
        f"Available columns: {list(table.columns)}"
    )


def optional_existing_column(
    table: pd.DataFrame,
    candidates: list[str],
) -> str | None:
    for column in candidates:
        if column in table.columns:
            return column

    return None


biology_summary = read_result_table(
    INPUT_FILES[
        "biology_summary"
    ]
)

biology_per_target = read_result_table(
    INPUT_FILES[
        "biology_per_target"
    ]
)

gradient_summary = read_result_table(
    INPUT_FILES[
        "gradient_summary"
    ]
)

gradient_per_target = read_result_table(
    INPUT_FILES[
        "gradient_per_target"
    ]
)

gradient_per_bin = read_result_table(
    INPUT_FILES[
        "gradient_per_bin"
    ]
)


optional_tables = {}

for key in [
    "biology_overall",
    "biology_paired",
    "biology_wins",
    "gradient_overall",
    "gradient_paired",
    "gradient_wins",
]:
    path = INPUT_FILES[
        key
    ]

    optional_tables[
        key
    ] = (
        read_result_table(
            path
        )
        if path.exists()
        else None
    )


BIO_MODEL = first_existing_column(
    biology_summary,
    [
        "model",
        "Model",
        "model_name",
    ],
    "biology_summary",
)

BIO_MAE = first_existing_column(
    biology_summary,
    [
        "MAE_mean",
        "mae_mean",
        "mean_MAE",
        "mean_mae",
        "MAE",
        "mae",
    ],
    "biology_summary",
)

BIO_MAE_SD = optional_existing_column(
    biology_summary,
    [
        "MAE_std",
        "mae_std",
        "MAE_sd",
        "mae_sd",
    ],
)

GRAD_MODEL = first_existing_column(
    gradient_summary,
    [
        "model",
        "Model",
        "model_name",
    ],
    "gradient_summary",
)

GRAD_MAE = first_existing_column(
    gradient_summary,
    [
        "MAE_mean",
        "mae_mean",
        "mean_MAE",
        "mean_mae",
        "MAE",
        "mae",
    ],
    "gradient_summary",
)

GRAD_MAE_SD = optional_existing_column(
    gradient_summary,
    [
        "MAE_std",
        "mae_std",
        "MAE_sd",
        "mae_sd",
    ],
)


BIO_PT_MODEL = first_existing_column(
    biology_per_target,
    [
        "model",
        "Model",
        "model_name",
    ],
    "biology_per_target",
)

BIO_PT_TARGET = first_existing_column(
    biology_per_target,
    [
        "gene_symbol",
        "perturbation",
        "target",
        "heldout_target",
    ],
    "biology_per_target",
)

BIO_PT_MAE = first_existing_column(
    biology_per_target,
    [
        "MAE",
        "mae",
        "target_MAE",
        "target_mae",
    ],
    "biology_per_target",
)

BIO_PT_PRED_MAG = optional_existing_column(
    biology_per_target,
    [
        "predicted_delta_magnitude",
        "predicted_magnitude",
        "pred_delta_magnitude",
    ],
)

BIO_PT_TRUE_MAG = optional_existing_column(
    biology_per_target,
    [
        "true_delta_magnitude",
        "true_magnitude",
        "observed_delta_magnitude",
    ],
)


GRAD_PT_MODEL = first_existing_column(
    gradient_per_target,
    [
        "model",
        "Model",
        "model_name",
    ],
    "gradient_per_target",
)

GRAD_PT_TARGET = first_existing_column(
    gradient_per_target,
    [
        "gene_symbol",
        "perturbation",
        "target",
        "heldout_target",
    ],
    "gradient_per_target",
)

GRAD_PT_MAE = first_existing_column(
    gradient_per_target,
    [
        "MAE",
        "mae",
        "target_MAE",
        "target_mae",
    ],
    "gradient_per_target",
)

GRAD_PT_PRED_MAG = first_existing_column(
    gradient_per_target,
    [
        "predicted_delta_magnitude",
        "predicted_magnitude",
        "pred_delta_magnitude",
    ],
    "gradient_per_target",
)

GRAD_PT_TRUE_MAG = first_existing_column(
    gradient_per_target,
    [
        "true_delta_magnitude",
        "true_magnitude",
        "observed_delta_magnitude",
    ],
    "gradient_per_target",
)


GRAD_BIN_MODEL = first_existing_column(
    gradient_per_bin,
    [
        "model",
        "Model",
        "model_name",
    ],
    "gradient_per_bin",
)

GRAD_BIN = first_existing_column(
    gradient_per_bin,
    [
        "kd_bin",
        "bin",
        "dose_bin",
        "strength_bin",
    ],
    "gradient_per_bin",
)

GRAD_BIN_PRED_MAG = first_existing_column(
    gradient_per_bin,
    [
        "predicted_delta_magnitude",
        "predicted_magnitude",
        "pred_delta_magnitude",
    ],
    "gradient_per_bin",
)

GRAD_BIN_TRUE_MAG = first_existing_column(
    gradient_per_bin,
    [
        "true_delta_magnitude",
        "true_magnitude",
        "observed_delta_magnitude",
    ],
    "gradient_per_bin",
)


print("\nBiology-informed summary:")
display(
    biology_summary
)

print("\nDose-aware summary:")
display(
    gradient_summary
)


# ============================================================
# 4. General plotting helpers
# ============================================================
def save_figure(
    figure: plt.Figure,
    filename: str,
) -> Path:
    path = (
        FIGURE_DIR
        / filename
    )

    figure.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    print(
        "Saved:",
        path,
    )

    return path


def baseline_name(
    models: list[str],
    preferred: list[str],
) -> str:
    for item in preferred:
        if item in models:
            return item

    return models[0]


def clean_model_label(
    label: str,
) -> str:
    replacements = {
        "Monotonic-TF-activity+weighted-loss": (
            "Monotonic TF activity\n+ weighted loss"
        ),
        "Monotonic-TF-activity": (
            "Monotonic TF activity"
        ),
        "Raw-TF-direct": (
            "Raw TF direct"
        ),
        "Shuffled-Hill-gradient": (
            "Shuffled Hill"
        ),
        "Linear-gradient": (
            "Linear gradient"
        ),
        "Strength-input": (
            "Strength input"
        ),
        "No-gradient": (
            "No gradient"
        ),
        "Hill-gradient": (
            "Hill gradient"
        ),
    }

    return replacements.get(
        str(
            label
        ),
        str(
            label
        ).replace(
            "-",
            " ",
        ),
    )


# ============================================================
# 5. Corrected sequential architecture
# ============================================================
fig, axis = plt.subplots(
    figsize=(
        16,
        6.5,
    )
)

axis.set_xlim(
    0,
    16,
)

axis.set_ylim(
    0,
    7,
)

axis.axis(
    "off"
)


def architecture_box(
    x: float,
    y: float,
    width: float,
    height: float,
    label: str,
    fontsize: int = 11,
) -> None:
    patch = FancyBboxPatch(
        (
            x,
            y,
        ),
        width,
        height,
        boxstyle=(
            "round,pad=0.04,"
            "rounding_size=0.12"
        ),
        fill=False,
        linewidth=1.6,
    )

    axis.add_patch(
        patch
    )

    axis.text(
        x + width / 2,
        y + height / 2,
        label,
        ha="center",
        va="center",
        fontsize=fontsize,
    )


def architecture_arrow(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
) -> None:
    axis.add_patch(
        FancyArrowPatch(
            (
                x1,
                y1,
            ),
            (
                x2,
                y2,
            ),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.5,
        )
    )


input_y = [
    5.2,
    3.8,
    2.4,
    1.0,
]

input_labels = [
    "Control-expression\nprofile",
    "ESM-2 target-protein\nembedding",
    "STRING physical-PPI\nfeatures",
    "Signed regulatory\nprior",
]

for y, label in zip(
    input_y,
    input_labels,
):
    architecture_box(
        0.4,
        y,
        2.6,
        0.9,
        label,
    )

architecture_box(
    4.3,
    2.65,
    2.6,
    1.4,
    "Feature\nintegration",
    fontsize=12,
)

architecture_box(
    8.0,
    2.65,
    2.5,
    1.4,
    "Residual network\npredicts maximum Δ",
    fontsize=12,
)

architecture_box(
    11.6,
    2.65,
    2.2,
    1.4,
    "Dose function\nscales Δ",
    fontsize=12,
)

architecture_box(
    14.6,
    2.65,
    1.1,
    1.4,
    "Final\nresponse",
    fontsize=11,
)

for y in input_y:
    architecture_arrow(
        3.0,
        y + 0.45,
        4.3,
        3.35,
    )

architecture_arrow(
    6.9,
    3.35,
    8.0,
    3.35,
)

architecture_arrow(
    10.5,
    3.35,
    11.6,
    3.35,
)

architecture_arrow(
    13.8,
    3.35,
    14.6,
    3.35,
)

axis.text(
    8.0,
    6.55,
    "Biology-informed, dose-aware virtual-cell model",
    ha="center",
    va="center",
    fontsize=19,
    fontweight="bold",
)

axis.text(
    12.7,
    1.7,
    "No-gradient: Δ is used directly\nLinear/Hill: Δ is dose-scaled",
    ha="center",
    va="center",
    fontsize=10,
)

axis.text(
    8.0,
    0.25,
    (
        "Perturbation identities are split before pseudobulk generation, "
        "preventing target leakage into validation or test data."
    ),
    ha="center",
    va="center",
    fontsize=10.5,
)

fig.tight_layout()

save_figure(
    fig,
    "model_architecture_polished.png",
)

plt.show()


# ============================================================
# 6. Clean MAE dot plot
# ============================================================
def plot_model_mae(
    summary: pd.DataFrame,
    model_col: str,
    mae_col: str,
    mae_sd_col: str | None,
    title: str,
    filename: str,
    baseline_preferences: list[str],
) -> tuple[pd.DataFrame, str]:
    plot_table = (
        summary[
            [
                model_col,
                mae_col,
            ]
            + (
                [
                    mae_sd_col,
                ]
                if mae_sd_col is not None
                else []
            )
        ]
        .dropna(
            subset=[
                model_col,
                mae_col,
            ]
        )
        .copy()
    )

    plot_table = (
        plot_table
        .groupby(
            model_col,
            as_index=False,
        )
        .agg({
            mae_col: "mean",
            **(
                {
                    mae_sd_col: "mean",
                }
                if mae_sd_col is not None
                else {}
            ),
        })
        .sort_values(
            mae_col,
            ascending=False,
        )
    )

    models = (
        plot_table[
            model_col
        ]
        .astype(str)
        .tolist()
    )

    baseline = baseline_name(
        models,
        baseline_preferences,
    )

    baseline_mae = float(
        plot_table.loc[
            plot_table[
                model_col
            ]
            .astype(str)
            .eq(
                baseline
            ),
            mae_col,
        ].iloc[0]
    )

    fig, axis = plt.subplots(
        figsize=(
            11.5,
            max(
                5.2,
                0.82
                * len(
                    plot_table
                ),
            ),
        )
    )

    y = np.arange(
        len(
            plot_table
        )
    )

    xerr = (
        plot_table[
            mae_sd_col
        ].to_numpy()
        if mae_sd_col is not None
        else None
    )

    axis.errorbar(
        plot_table[
            mae_col
        ],
        y,
        xerr=xerr,
        fmt="o",
        capsize=4,
        markersize=8,
        zorder=3,
    )

    axis.axvline(
        baseline_mae,
        linestyle="--",
        linewidth=1.2,
        label=(
            f"{clean_model_label(baseline)} reference"
        ),
    )

    axis.set_yticks(
        y
    )

    axis.set_yticklabels([
        clean_model_label(
            value
        )
        for value in plot_table[
            model_col
        ]
    ])

    axis.set_xlabel(
        "Mean absolute error"
    )

    axis.set_ylabel(
        "Model"
    )

    axis.set_title(
        title,
        pad=14,
    )

    axis.legend(
        frameon=False,
        loc="lower left",
    )

    axis.grid(
        axis="x",
        alpha=0.25,
    )

    values = plot_table[
        mae_col
    ].to_numpy(
        dtype=float
    )

    if xerr is not None:
        left_extent = float(
            np.min(
                values
                - xerr
            )
        )
        right_extent = float(
            np.max(
                values
                + xerr
            )
        )
    else:
        left_extent = float(
            values.min()
        )
        right_extent = float(
            values.max()
        )

    span = max(
        right_extent
        - left_extent,
        1e-4,
    )

    text_x = (
        right_extent
        + 0.06
        * span
    )

    axis.set_xlim(
        left_extent
        - 0.05
        * span,
        right_extent
        + 0.42
        * span,
    )

    axis.axvline(
        text_x
        - 0.025
        * span,
        linewidth=0.8,
        alpha=0.25,
    )

    axis.text(
        text_x,
        len(
            plot_table
        )
        - 0.35,
        "MAE   Δ vs baseline",
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )

    for index, row in enumerate(
        plot_table.itertuples(
            index=False
        )
    ):
        value = float(
            getattr(
                row,
                mae_col,
            )
        )

        improvement = (
            baseline_mae
            - value
        )

        axis.text(
            text_x,
            index,
            (
                f"{value:.4f}   "
                f"{improvement:+.4f}"
            ),
            ha="left",
            va="center",
            fontsize=9,
        )

    fig.tight_layout()

    save_figure(
        fig,
        filename,
    )

    plt.show()

    return (
        plot_table,
        baseline,
    )


biology_plot, biology_baseline = plot_model_mae(
    summary=biology_summary,
    model_col=BIO_MODEL,
    mae_col=BIO_MAE,
    mae_sd_col=BIO_MAE_SD,
    title=(
        "Biology-informed models on held-out perturbations"
    ),
    filename=(
        "biology_prior_model_comparison_polished.png"
    ),
    baseline_preferences=[
        "Baseline",
        "ESM-only",
        "No-change",
    ],
)


gradient_plot, gradient_baseline = plot_model_mae(
    summary=gradient_summary,
    model_col=GRAD_MODEL,
    mae_col=GRAD_MAE,
    mae_sd_col=GRAD_MAE_SD,
    title=(
        "Dose-aware models on held-out perturbations"
    ),
    filename=(
        "dose_model_comparison_polished.png"
    ),
    baseline_preferences=[
        "No-gradient",
        "Baseline",
        "No-change",
    ],
)


# ============================================================
# 7. Ranked held-out target performance
# ============================================================
available_gradient_models = (
    gradient_per_target[
        GRAD_PT_MODEL
    ]
    .dropna()
    .astype(str)
    .unique()
    .tolist()
)

preferred_gradient_model = (
    "Hill-gradient"
    if "Hill-gradient"
    in available_gradient_models
    else available_gradient_models[0]
)


target_table = (
    gradient_per_target.loc[
        gradient_per_target[
            GRAD_PT_MODEL
        ]
        .astype(str)
        .eq(
            preferred_gradient_model
        )
    ]
    .groupby(
        GRAD_PT_TARGET,
        as_index=False,
    )
    .agg(
        mean_MAE=(
            GRAD_PT_MAE,
            "mean",
        ),
    )
    .sort_values(
        "mean_MAE",
        ascending=True,
    )
)


fig, axis = plt.subplots(
    figsize=(
        10,
        max(
            7,
            0.33
            * len(
                target_table
            ),
        ),
    )
)

y = np.arange(
    len(
        target_table
    )
)

axis.hlines(
    y,
    xmin=0,
    xmax=target_table[
        "mean_MAE"
    ],
    linewidth=1.5,
)

axis.plot(
    target_table[
        "mean_MAE"
    ],
    y,
    "o",
)

axis.set_yticks(
    y
)

axis.set_yticklabels(
    target_table[
        GRAD_PT_TARGET
    ]
)

axis.set_xlabel(
    "Mean absolute error"
)

axis.set_ylabel(
    "Held-out target"
)

axis.set_title(
    (
        "Generalization across unseen perturbations — "
        f"{clean_model_label(preferred_gradient_model)}"
    )
)

axis.grid(
    axis="x",
    alpha=0.25,
)

fig.tight_layout()

save_figure(
    fig,
    "heldout_target_performance_polished.png",
)

plt.show()


# ============================================================
# 8. Predicted versus observed magnitude across dose models
# ============================================================
calibration_table = (
    gradient_per_target[
        [
            GRAD_PT_MODEL,
            GRAD_PT_TARGET,
            GRAD_PT_PRED_MAG,
            GRAD_PT_TRUE_MAG,
        ]
    ]
    .dropna()
    .groupby(
        [
            GRAD_PT_MODEL,
            GRAD_PT_TARGET,
        ],
        as_index=False,
    )
    .agg(
        predicted_magnitude=(
            GRAD_PT_PRED_MAG,
            "mean",
        ),
        observed_magnitude=(
            GRAD_PT_TRUE_MAG,
            "mean",
        ),
    )
)


fig, axis = plt.subplots(
    figsize=(
        8.5,
        7.5,
    )
)

for model, group in calibration_table.groupby(
    GRAD_PT_MODEL
):
    axis.scatter(
        group[
            "observed_magnitude"
        ],
        group[
            "predicted_magnitude"
        ],
        alpha=0.55,
        s=45,
        label=clean_model_label(
            model
        ),
    )


maximum = float(
    max(
        calibration_table[
            "observed_magnitude"
        ].max(),
        calibration_table[
            "predicted_magnitude"
        ].max(),
    )
)

axis.plot(
    [
        0,
        maximum,
    ],
    [
        0,
        maximum,
    ],
    linestyle="--",
    linewidth=1.5,
    label="Perfect agreement",
)

axis.set_xlim(
    0,
    maximum
    * 1.05,
)

axis.set_ylim(
    0,
    maximum
    * 1.05,
)

axis.set_xlabel(
    "Observed response magnitude"
)

axis.set_ylabel(
    "Predicted response magnitude"
)

axis.set_title(
    "Response-magnitude calibration across held-out targets"
)

axis.legend(
    frameon=False,
    bbox_to_anchor=(
        1.02,
        1,
    ),
    loc="upper left",
)

axis.grid(
    alpha=0.25,
)

fig.tight_layout()

save_figure(
    fig,
    "predicted_vs_observed_magnitude_by_model.png",
)

plt.show()


# ============================================================
# 9. Calibration error and correlation by model
# ============================================================
calibration_rows = []

for model, group in calibration_table.groupby(
    GRAD_PT_MODEL
):
    observed = group[
        "observed_magnitude"
    ].to_numpy(
        dtype=float
    )

    predicted = group[
        "predicted_magnitude"
    ].to_numpy(
        dtype=float
    )

    magnitude_mae = float(
        np.mean(
            np.abs(
                predicted
                - observed
            )
        )
    )

    if (
        np.std(
            observed
        ) > 0
        and np.std(
            predicted
        ) > 0
    ):
        correlation = float(
            np.corrcoef(
                observed,
                predicted,
            )[
                0,
                1,
            ]
        )
    else:
        correlation = np.nan

    calibration_rows.append({
        "model": model,
        "magnitude_MAE": (
            magnitude_mae
        ),
        "magnitude_correlation": (
            correlation
        ),
        "predicted_to_observed_ratio": float(
            predicted.mean()
            / (
                observed.mean()
                + 1e-12
            )
        ),
    })


calibration_summary = (
    pd.DataFrame(
        calibration_rows
    )
    .sort_values(
        "magnitude_MAE",
        ascending=False,
    )
)


fig, axis = plt.subplots(
    figsize=(
        9,
        max(
            4.8,
            0.65
            * len(
                calibration_summary
            ),
        ),
    )
)

y = np.arange(
    len(
        calibration_summary
    )
)

axis.barh(
    y,
    calibration_summary[
        "magnitude_MAE"
    ],
)

axis.set_yticks(
    y
)

axis.set_yticklabels([
    clean_model_label(
        value
    )
    for value in calibration_summary[
        "model"
    ]
])

axis.set_xlabel(
    "Mean |predicted − observed| response magnitude"
)

axis.set_ylabel(
    "Model"
)

axis.set_title(
    "Response-magnitude calibration error"
)

axis.grid(
    axis="x",
    alpha=0.25,
)

for index, row in calibration_summary.reset_index(
    drop=True
).iterrows():
    axis.text(
        row[
            "magnitude_MAE"
        ],
        index,
        (
            f" {row['magnitude_MAE']:.4f}; "
            f"r={row['magnitude_correlation']:.2f}; "
            f"ratio={row['predicted_to_observed_ratio']:.2f}"
        ),
        va="center",
        fontsize=9,
    )

fig.tight_layout()

save_figure(
    fig,
    "magnitude_calibration_error_by_model.png",
)

plt.show()


calibration_summary.to_csv(
    RESULT_DIR
    / "magnitude_calibration_summary.csv",
    index=False,
)


# ============================================================
# 10. Observed and predicted dose trajectories among models
# ============================================================
trajectory_table = (
    gradient_per_bin[
        [
            GRAD_BIN_MODEL,
            GRAD_BIN,
            GRAD_BIN_PRED_MAG,
            GRAD_BIN_TRUE_MAG,
        ]
    ]
    .dropna()
    .groupby(
        [
            GRAD_BIN_MODEL,
            GRAD_BIN,
        ],
        as_index=False,
    )
    .agg(
        predicted_magnitude=(
            GRAD_BIN_PRED_MAG,
            "mean",
        ),
        observed_magnitude=(
            GRAD_BIN_TRUE_MAG,
            "mean",
        ),
    )
)


observed_trajectory = (
    trajectory_table
    .groupby(
        GRAD_BIN,
        as_index=False,
    )
    .agg(
        observed_magnitude=(
            "observed_magnitude",
            "mean",
        )
    )
    .sort_values(
        GRAD_BIN
    )
)


fig, axis = plt.subplots(
    figsize=(
        9,
        6,
    )
)

axis.plot(
    observed_trajectory[
        GRAD_BIN
    ],
    observed_trajectory[
        "observed_magnitude"
    ],
    marker="o",
    linewidth=3,
    label="Observed",
)

for model, group in trajectory_table.groupby(
    GRAD_BIN_MODEL
):
    group = group.sort_values(
        GRAD_BIN
    )

    axis.plot(
        group[
            GRAD_BIN
        ],
        group[
            "predicted_magnitude"
        ],
        marker="o",
        linewidth=1.7,
        label=(
            f"Predicted — "
            f"{clean_model_label(model)}"
        ),
    )

axis.set_xlabel(
    "Knockdown-strength bin"
)

axis.set_ylabel(
    "Mean absolute transcriptomic response"
)

axis.set_title(
    "Observed and predicted dose-response trajectories"
)

axis.set_xticks(
    sorted(
        trajectory_table[
            GRAD_BIN
        ].unique()
    )
)

axis.legend(
    frameon=False,
    bbox_to_anchor=(
        1.02,
        1,
    ),
    loc="upper left",
)

axis.grid(
    alpha=0.25,
)

fig.tight_layout()

save_figure(
    fig,
    "dose_response_trajectories_all_models.png",
)

plt.show()


# ============================================================
# 11. Optional biology-model calibration
# ============================================================
if (
    BIO_PT_PRED_MAG is not None
    and BIO_PT_TRUE_MAG is not None
):
    biology_calibration = (
        biology_per_target[
            [
                BIO_PT_MODEL,
                BIO_PT_TARGET,
                BIO_PT_PRED_MAG,
                BIO_PT_TRUE_MAG,
            ]
        ]
        .dropna()
        .groupby(
            [
                BIO_PT_MODEL,
                BIO_PT_TARGET,
            ],
            as_index=False,
        )
        .agg(
            predicted_magnitude=(
                BIO_PT_PRED_MAG,
                "mean",
            ),
            observed_magnitude=(
                BIO_PT_TRUE_MAG,
                "mean",
            ),
        )
    )

    fig, axis = plt.subplots(
        figsize=(
            8.5,
            7.5,
        )
    )

    for model, group in biology_calibration.groupby(
        BIO_PT_MODEL
    ):
        axis.scatter(
            group[
                "observed_magnitude"
            ],
            group[
                "predicted_magnitude"
            ],
            alpha=0.55,
            s=45,
            label=clean_model_label(
                model
            ),
        )

    maximum = float(
        max(
            biology_calibration[
                "observed_magnitude"
            ].max(),
            biology_calibration[
                "predicted_magnitude"
            ].max(),
        )
    )

    axis.plot(
        [
            0,
            maximum,
        ],
        [
            0,
            maximum,
        ],
        linestyle="--",
        linewidth=1.5,
        label="Perfect agreement",
    )

    axis.set_xlim(
        0,
        maximum
        * 1.05,
    )

    axis.set_ylim(
        0,
        maximum
        * 1.05,
    )

    axis.set_xlabel(
        "Observed response magnitude"
    )

    axis.set_ylabel(
        "Predicted response magnitude"
    )

    axis.set_title(
        "Biology-model response-magnitude calibration"
    )

    axis.legend(
        frameon=False,
        bbox_to_anchor=(
            1.02,
            1,
        ),
        loc="upper left",
    )

    axis.grid(
        alpha=0.25,
    )

    fig.tight_layout()

    save_figure(
        fig,
        "biology_models_predicted_vs_observed_magnitude.png",
    )

    plt.show()


# ============================================================
# 12. Clean README overview infographic
# ============================================================
best_biology_row = (
    biology_summary
    .sort_values(
        BIO_MAE,
        ascending=True,
    )
    .iloc[0]
)

best_gradient_row = (
    gradient_summary
    .sort_values(
        GRAD_MAE,
        ascending=True,
    )
    .iloc[0]
)

best_calibration_row = (
    calibration_summary
    .sort_values(
        "magnitude_MAE",
        ascending=True,
    )
    .iloc[0]
)


fig, axis = plt.subplots(
    figsize=(
        17,
        10.5,
    )
)

axis.set_xlim(
    0,
    17,
)

axis.set_ylim(
    0,
    10.5,
)

axis.axis(
    "off"
)


def wrap_multiline_text(
    text: str,
    width: int,
) -> str:
    wrapped_lines = []

    for line in str(
        text
    ).splitlines():
        stripped = line.strip()

        if not stripped:
            wrapped_lines.append(
                ""
            )
            continue

        bullet_prefix = ""

        if stripped.startswith(
            "• "
        ):
            bullet_prefix = "• "
            stripped = stripped[2:]

        wrapped = textwrap.wrap(
            stripped,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )

        if not wrapped:
            wrapped_lines.append(
                bullet_prefix
            )
            continue

        wrapped_lines.append(
            bullet_prefix
            + wrapped[0]
        )

        for continuation in wrapped[1:]:
            wrapped_lines.append(
                (
                    "  "
                    if bullet_prefix
                    else ""
                )
                + continuation
            )

    return "\n".join(
        wrapped_lines
    )


def info_box(
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    wrap_width: int = 44,
    body_fontsize: float = 9.3,
) -> None:
    patch = FancyBboxPatch(
        (
            x,
            y,
        ),
        width,
        height,
        boxstyle=(
            "round,pad=0.05,"
            "rounding_size=0.12"
        ),
        linewidth=1.4,
        fill=False,
    )

    axis.add_patch(
        patch
    )

    axis.text(
        x + 0.28,
        y + height - 0.28,
        title,
        ha="left",
        va="top",
        fontsize=12.2,
        fontweight="bold",
    )

    axis.text(
        x + 0.28,
        y + height - 0.86,
        wrap_multiline_text(
            body,
            width=wrap_width,
        ),
        ha="left",
        va="top",
        fontsize=body_fontsize,
        linespacing=1.28,
    )


axis.text(
    8.5,
    10.05,
    "Predicting transcriptional responses to unseen perturbations",
    ha="center",
    va="center",
    fontsize=21,
    fontweight="bold",
)

axis.text(
    8.5,
    9.55,
    (
        "A leakage-controlled virtual-cell workflow integrating "
        "protein sequence, network priors, regulatory information, "
        "and perturbation strength."
    ),
    ha="center",
    va="center",
    fontsize=11,
)

info_box(
    0.55,
    5.55,
    5.1,
    3.25,
    "Problem and evaluation",
    (
        "• Predict expression responses for perturbation targets "
        "absent from training\n"
        "• Split perturbation identities before pseudobulk generation\n"
        "• Evaluate MAE, directional similarity, rank correlation, "
        "top-gene overlap, and response magnitude"
    ),
    wrap_width=46,
)

info_box(
    5.95,
    5.55,
    5.1,
    3.25,
    "Model design",
    (
        "• Control-expression profile\n"
        "• ESM-2 target-protein embedding\n"
        "• STRING physical-PPI features\n"
        "• Signed regulatory prior\n"
        "• Residual neural network\n"
        "• Linear and Hill dose transformations"
    ),
    wrap_width=43,
)

info_box(
    11.35,
    5.55,
    5.1,
    3.25,
    "Controls and analysis",
    (
        "• Baseline and no-gradient comparisons\n"
        "• Shuffled-Hill control\n"
        "• Multiple held-out perturbation targets\n"
        "• Per-target error analysis\n"
        "• Predicted-versus-observed calibration\n"
        "• Dose-bin trajectory comparison"
    ),
    wrap_width=43,
)

info_box(
    0.55,
    1.25,
    5.1,
    3.45,
    "Best biology-informed MAE",
    (
        f"Model: {clean_model_label(best_biology_row[BIO_MODEL])}\n\n"
        f"Mean MAE: {float(best_biology_row[BIO_MAE]):.4f}\n\n"
        "This result belongs to the biology-prior experiment and "
        "should be interpreted against its own baseline."
    ),
    wrap_width=42,
)

info_box(
    5.95,
    1.25,
    5.1,
    3.45,
    "Best dose-aware MAE",
    (
        f"Model: {clean_model_label(best_gradient_row[GRAD_MODEL])}\n\n"
        f"Mean MAE: {float(best_gradient_row[GRAD_MAE]):.4f}\n\n"
        "Compare this value with the no-gradient reference to assess "
        "the contribution of perturbation strength."
    ),
    wrap_width=42,
)

info_box(
    11.35,
    1.25,
    5.1,
    3.45,
    "Best magnitude calibration",
    (
        f"Model: {clean_model_label(best_calibration_row['model'])}\n\n"
        f"Magnitude MAE: "
        f"{float(best_calibration_row['magnitude_MAE']):.4f}\n"
        f"Correlation: "
        f"{float(best_calibration_row['magnitude_correlation']):.2f}\n"
        f"Predicted/observed ratio: "
        f"{float(best_calibration_row['predicted_to_observed_ratio']):.2f}"
    ),
    wrap_width=40,
)

axis.text(
    8.5,
    0.45,
    (
        "Biology-informed and dose-aware experiments are summarized "
        "separately because their raw MAE values are not automatically "
        "one shared leaderboard."
    ),
    ha="center",
    va="center",
    fontsize=9.5,
)

fig.tight_layout()

save_figure(
    fig,
    "portfolio_summary_polished.png",
)

plt.show()


# ============================================================
# 13. Export GitHub tables and README section
# ============================================================
biology_export = biology_summary.copy()

biology_export.insert(
    0,
    "experiment",
    "biology_informed",
)

gradient_export = gradient_summary.copy()

gradient_export.insert(
    0,
    "experiment",
    "dose_aware",
)

all_columns = sorted(
    set(
        biology_export.columns
    )
    | set(
        gradient_export.columns
    )
)

combined_summary = pd.concat(
    [
        biology_export.reindex(
            columns=all_columns
        ),
        gradient_export.reindex(
            columns=all_columns
        ),
    ],
    ignore_index=True,
)

combined_summary.to_csv(
    RESULT_DIR
    / "model_summary.csv",
    index=False,
)


readme_text = f"""## Model overview

![Model architecture](results/figures/model_architecture_polished.png)

The workflow predicts transcriptional responses to perturbations that
were excluded from training. Perturbation identities are split before
pseudobulk construction to prevent target leakage.

## Project summary

![Portfolio summary](results/figures/portfolio_summary_polished.png)

## Biology-informed model comparison

![Biology-informed comparison](results/figures/biology_prior_model_comparison_polished.png)

The lowest mean MAE in this experiment was obtained by
**{clean_model_label(best_biology_row[BIO_MODEL])}**
(**{float(best_biology_row[BIO_MAE]):.4f}**).

## Dose-aware model comparison

![Dose-aware comparison](results/figures/dose_model_comparison_polished.png)

The lowest mean MAE in this experiment was obtained by
**{clean_model_label(best_gradient_row[GRAD_MODEL])}**
(**{float(best_gradient_row[GRAD_MAE]):.4f}**).

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
"""


(
    RESULT_DIR
    / "README_results_section.md"
).write_text(
    readme_text,
    encoding="utf-8",
)


print("\nGenerated polished figures:")

for path in sorted(
    FIGURE_DIR.glob(
        "*polished.png"
    )
):
    print(
        " -",
        path.name,
    )

for path in [
    FIGURE_DIR
    / "predicted_vs_observed_magnitude_by_model.png",
    FIGURE_DIR
    / "magnitude_calibration_error_by_model.png",
    FIGURE_DIR
    / "dose_response_trajectories_all_models.png",
]:
    if path.exists():
        print(
            " -",
            path.name,
        )


print(
    "\nSaved README section:",
    RESULT_DIR
    / "README_results_section.md",
)
