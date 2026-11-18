"""Recover paired left-thalamus fields and audit the equal-cortical-P99 solve."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import simnibs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from analyze import weighted_quantile
from define_target import topology_hash

ROOT = HERE.parents[1]
COHORT = ("sub01", "sub09", "sub10", "sub12", "sub15")
PILOT = "FieldTrip_Subject01"
COILS = ("h1", "double_cone", "figure8")
P99 = 30.0
DRIVE = 1e8
OUT = HERE / "results"


def save_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=tuple(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def load_csv(path):
    with path.open(encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def analyze_subject(subject):
    pilot = subject == PILOT
    cache = ROOT / "tmp" / ("brain_field" if pilot else "brain_cohort")
    cache = cache if pilot else cache / subject
    records = HERE.parent / "results" if pilot else OUT / subject
    target = json.loads((records / "target.json").read_text(encoding="utf-8"))
    solver = json.loads((records / "solver_run.json").read_text(encoding="utf-8"))
    pose = json.loads((records / "placement.json").read_text(encoding="utf-8"))
    assert solver["didt_A_per_s"] == DRIVE and target["radius_mm"] == 5
    file = records / ("summary.csv" if pilot else "matched_check.csv")
    checks = {r["coil"]: r for r in load_csv(file)}
    assert set(checks) == set(COILS)
    rows, radii = [], []
    with np.load(cache / "target_roi.npz") as roi:
        ids = roi["element_numbers"] - 1
        weights = roi["volumes_mm3"]
        cortex_ids = roi["reference_element_numbers"] - 1
        cortex_weights = roi["reference_volumes_mm3"]
        regional = roi["regional_element_numbers"] - 1
        region_weights = roi["regional_volumes_mm3"]
        distances = roi["distances_mm"]
        assert len(ids) == target["roi_element_count"] and np.all(weights > 0)
        assert len(cortex_ids) == target["cortical_element_count"] and np.all(cortex_weights > 0)
        for coil in COILS:
            mesh = simnibs.read_msh(str(cache / "fields" / f"{coil}.msh"))
            if topology_hash(mesh) != target["topology_sha256"]:
                raise ValueError(f"{subject}: field and anatomical mesh do not match")
            centers = mesh.nodes.node_coord[mesh.elm.node_number_list[ids] - 1].mean(axis=1)
            np.testing.assert_allclose(centers, roi["centroids_mm"], rtol=0, atol=1e-5)
            E = mesh.field["E"].value
            brain = np.linalg.norm(E[ids], axis=1)
            cortex = np.linalg.norm(E[cortex_ids], axis=1)
            if not np.all(np.isfinite(cortex)) or not np.all(np.isfinite(brain)):
                raise ValueError("Unexpected nonfinite brain field")
            fixed = np.average(brain, weights=weights)
            p99 = weighted_quantile(cortex, cortex_weights, .99)
            if p99 <= 0:
                raise ValueError("Invalid reference P99")
            gain = P99 / p99
            expected = np.average(np.linalg.norm(gain * E[ids], axis=1), weights=weights)
            scaled_p99 = weighted_quantile(gain * cortex, cortex_weights, .99)
            np.testing.assert_allclose(scaled_p99, P99, atol=1e-10, rtol=0)
            check = checks[coil]
            if pilot:
                np.testing.assert_allclose([fixed, p99],
                    [float(check["E_roi_mean_V_m"]), float(check["E_reference_p99_V_m"])], rtol=1e-10)
                matched, matched_ref = expected, scaled_p99
            else:
                matched, matched_ref = float(check["matched_roi_mean_V_m"]), float(check["matched_p99_V_m"])
                np.testing.assert_allclose([fixed, p99, gain, DRIVE*gain, expected, P99],
                    [float(check["fixed_roi_mean_V_m"]), float(check["fixed_cortex_p99_V_m"]),
                     float(check["gain"]), float(check["matched_didt_A_per_s"]), matched, matched_ref],
                    rtol=1e-6, atol=1e-8)
            group = "pilot" if pilot else "replication"
            rows.append(dict(subject=subject, group=group, coil=coil, ROI_radius_mm=5,
                ROI_elements=len(ids), ROI_volume_mm3=weights.sum(), upper_scalp_depth_mm=pose["superior_depth_mm"],
                wire_gap_mm=pose["coils"][coil]["wire_gap_mm"], fixed_didt_A_per_us=DRIVE/1e6,
                fixed_target_V_m=fixed, fixed_cortex_P99_V_m=p99, gain=gain,
                matched_didt_A_per_us=DRIVE*gain/1e6, matched_cortex_P99_V_m=matched_ref,
                matched_target_V_m=matched, matched_method=("field_scaling" if pilot else "second_FEM_solve")))
            regional_norm = np.linalg.norm(E[regional], axis=1)
            for radius in (3.0, 5.0, 7.5):
                mask = distances <= radius
                mean = np.average(regional_norm[mask], weights=region_weights[mask])
                radii.append(dict(subject=subject, group=group, coil=coil, radius_mm=radius,
                    elements=int(mask.sum()), volume_mm3=region_weights[mask].sum(),
                    fixed_target_V_m=mean, scaled_to_P99_30_target_V_m=mean*gain))
            del mesh, E
    save_csv(OUT / subject / "comparison.csv", rows)
    save_csv(OUT / subject / "radius_comparison.csv", radii)


def aggregate():
    rows, radii, paired = [], [], []
    for subject in (*COHORT, PILOT):
        data = load_csv(OUT / subject / "comparison.csv")
        actual = {r["coil"]: r for r in data}
        assert set(actual) == set(COILS) and len(data) == len(COILS)
        rows.extend(data)
        radii.extend(load_csv(OUT / subject / "radius_comparison.csv"))
        for drive in ("fixed", "matched"):
            scores = {c: float(actual[c][drive + "_target_V_m"]) for c in COILS}
            paired.append(dict(subject=subject, group=data[0]["group"], condition=drive,
                               h1_over_double_cone=scores["h1"]/scores["double_cone"],
                               h1_over_figure8=scores["h1"]/scores["figure8"],
                               h1_minus_double_cone_V_m=scores["h1"]-scores["double_cone"],
                               h1_minus_figure8_V_m=scores["h1"]-scores["figure8"],
                               order=" > ".join(sorted(scores, key=scores.get, reverse=True))))
    summaries = []
    for coil in COILS:
        for condition in ("fixed", "matched"):
            scores = [float(r[condition + "_target_V_m"]) for r in rows if r["group"] == "replication" and r["coil"] == coil]
            summaries.append(dict(coil=coil, condition=condition, n_individuals=len(scores),
                min_target_V_m=min(scores), median_target_V_m=np.median(scores), max_target_V_m=max(scores)))
    save_csv(OUT / "per_subject.csv", rows)
    save_csv(OUT / "paired_comparisons.csv", paired)
    save_csv(OUT / "radius_sensitivity.csv", radii)
    save_csv(OUT / "cohort_summary.csv", summaries)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject", choices=(*COHORT, PILOT, "aggregate"))
    subject = parser.parse_args().subject
    if subject == "aggregate":
        aggregate()
    else:
        analyze_subject(subject)
