"""Summarize the five subjects and compare the re-solved target condition with the scaled-field volumes."""
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "results"
IDS = {"sub01", "sub09", "sub10", "sub12", "sub15"}
COILS = {"h1", "double_cone", "figure8"}
CONDITIONS = {"cortex_P99_30", "target_mean_1"}
KEY = ["subject", "coil", "condition"]


def main():
    conditions = pd.read_csv(OUT / "conditions.csv")
    thresholds = pd.read_csv(OUT / "thresholds.csv")
    assert set(conditions.subject) == IDS and set(conditions.coil) == COILS
    assert set(conditions.condition) == CONDITIONS and len(conditions) == 30
    assert not conditions.duplicated(KEY).any()
    assert not thresholds.duplicated([*KEY, "threshold_V_m"]).any()
    assert len(thresholds) == 210
    threshold_metrics = ["cortex_cm3", "cortex_percent", "cortex_elements", "brain_cm3",
                         "brain_percent", "off_target_cm3", "off_target_percent"]
    condition_metrics = ["didt_A_per_us", "cortex_P99_V_m", "target_mean_V_m",
                         "cortex_half_P99_cm3", "cortex_half_P99_percent"]
    volumes = thresholds.melt([*KEY, "threshold_V_m"], threshold_metrics, var_name="metric")
    readouts = conditions.melt(KEY, condition_metrics, var_name="metric")
    rows = pd.concat([volumes, readouts], ignore_index=True)
    summary = rows.groupby(["condition", "threshold_V_m", "metric", "coil"], dropna=False).value.agg(
        n_subjects="size", minimum="min", median="median", maximum="max").reset_index()
    assert summary.n_subjects.eq(len(IDS)).all()

    direct = pd.concat([pd.read_csv(OUT / "verification" / f"{s}.csv") for s in sorted(IDS)])
    condition = conditions[conditions.condition == "target_mean_1"]
    check = direct.merge(condition, on=["subject", "coil"], validate="one_to_one", suffixes=("_solve", "_scaled"))
    assert len(check) == len(IDS) * len(COILS)
    for metric in ("didt_A_per_us", "cortex_P99_V_m", "target_mean_V_m"):
        np.testing.assert_allclose(check[metric+"_solve"], check[metric+"_scaled"], rtol=1e-8, atol=1e-10)
    for threshold in [1, 15]:
        selected = thresholds[(thresholds.condition == "target_mean_1") & (thresholds.threshold_V_m == threshold)]
        comparison = direct.merge(selected, on=["subject", "coil"], validate="one_to_one")
        np.testing.assert_allclose(comparison[f"cortex_ge_{threshold}_cm3"], comparison.cortex_cm3, atol=1e-8, rtol=0)
        assert comparison[f"changed_elements_at_{threshold}"].eq(0).all()
    assert direct.ksp_converged_reason.gt(0).all()
    summary.to_csv(OUT / "cohort_summary.csv", index=False)
    print(f"{len(direct)} additional FEM solves. Max |target-1| = "
          f"{(direct.target_mean_V_m-1).abs().max():.3e} V/m; "
          f"max cortical vector relative L2 error = {direct.cortex_vector_relative_L2.max():.3e}.")


if __name__ == "__main__":
    main()
