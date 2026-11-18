"""Compute tissue-volume exceedance curves on five paired native TMS head meshes.

Run from the project root with simnibs_python brain_field/cohort/focality/measure.py.
All inputs are the fixed-drive FEM vector fields and anatomies of the preceding
cohort experiment. The two prescribed reference levels specify positive linear
scalings of each subject-coil field.
"""
import hashlib
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
import simnibs

HERE = Path(__file__).resolve().parent
COHORT = HERE.parent
sys.path.insert(0, str(COHORT))
import measure_cohort as prior
from define_target import topology_hash
from analyze import weighted_quantile

CACHE = prior.ROOT / "tmp" / "brain_cohort"
RESULTS = HERE / "results"
THRESHOLDS = np.array([1.0, 7.5, 15.0, 22.5, 30.0, 45.0, 60.0])
CONDITIONS = {"cortex_P99_30": 30.0, "target_mean_1": 1.0}


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def coordinate_hash(mesh):
    return hashlib.sha256(np.ascontiguousarray(mesh.nodes.node_coord, dtype="<f8")).hexdigest()


def tetra_volumes(mesh, element_ids):
    """Native mm coordinates give cm³ after dividing the scalar triple by 6000."""
    volumes = np.empty(len(element_ids))
    for start in range(0, len(element_ids), 60000):
        ids = element_ids[start:start + 60000]
        vertices = mesh.nodes.node_coord[mesh.elm.node_number_list[ids, :4] - 1]
        a, b, c = (vertices[:, j] - vertices[:, 0] for j in (1, 2, 3))
        volumes[start:start + len(ids)] = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6000
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0):
        raise ValueError("Nonpositive or invalid native tissue volume")
    return volumes


class Exceedance:
    def __init__(self, values, volumes):
        order = np.argsort(values)
        self.fields = values[order]
        # Reverse accumulation avoids subtracting nearly identical full-volume sums.
        self.sums = np.r_[np.cumsum(volumes[order][::-1], dtype=np.float64)[::-1], 0.0]

    def above(self, cutoff):
        return self.sums[np.searchsorted(self.fields, cutoff, side="left")]

    def element_counts(self, cutoff):
        return len(self.fields) - np.searchsorted(self.fields, cutoff, side="left")


def volume_rows(subject, coil, condition, thresholds, source_reference, reference_level, curves, totals):
    # An exact reference-level threshold maps to its own source element, including
    # the tetrahedron that establishes the discrete weighted P99.
    cutoff = thresholds / reference_level * source_reference
    cortex = curves["cortex"].above(cutoff)
    brain = curves["brain"].above(cutoff)
    target = curves["target"].above(cutoff)
    cortical_count = curves["cortex"].element_counts(cutoff)
    brain_count = curves["brain"].element_counts(cutoff)
    target_count = curves["target"].element_counts(cutoff)
    off_target = brain - target
    for values in (cortex, brain, off_target):
        if np.any(np.diff(values) > 1e-9) or np.min(values) < -1e-10:
            raise ValueError("Invalid threshold-volume curve")
    if np.any(cortex > off_target + 1e-8):
        raise ValueError("Cortical comparator must be inside the off-target brain mask")
    return [
        dict(subject=subject, coil=coil, condition=condition, threshold_V_m=float(t),
             cortex_cm3=c, cortex_percent=100*c/totals["cortex"],
             brain_cm3=b, brain_percent=100*b/totals["brain"],
             off_target_cm3=o, off_target_percent=100*o/totals["off_target"],
             cortex_elements=nc, brain_elements=nb, off_target_elements=nb-nt)
        for t, c, b, o, nc, nb, nt in zip(thresholds, cortex, brain, off_target,
                                           cortical_count, brain_count, target_count)
    ]


def measure_subject(subject, checks, manifest):
    folder = COHORT / "results" / subject
    target = json.loads((folder / "target.json").read_text())
    solve = json.loads((folder / "solver_run.json").read_text())
    roi_file = CACHE / subject / "target_roi.npz"
    head_file = CACHE / "heads" / subject / f"{subject}.msh"
    assert solve["subject"] == target["subject"] == subject and target["radius_mm"] == 5
    assert solve["didt_A_per_s"] == prior.DRIVE
    assert sha256(roi_file) == solve["roi_file_sha256"]
    assert sha256(head_file) == solve["mesh_sha256"] == target["mesh_sha256"] == manifest[head_file.name]["sha256"]
    assert sha256(folder / "placement.json") == solve["placement_sha256"]
    with np.load(roi_file) as roi:
        cortex_ids = roi["reference_element_numbers"] - 1
        target_ids = roi["element_numbers"] - 1
        cortex_weights = roi["reference_volumes_mm3"]
        target_weights = roi["volumes_mm3"]
        target_centers = roi["centroids_mm"]
    head = simnibs.read_msh(str(head_file))
    node_digest = coordinate_hash(head)
    assert topology_hash(head) == target["topology_sha256"]
    for ids in (cortex_ids, target_ids):
        assert np.all(np.diff(ids) > 0)
        assert np.all((head.elm.elm_type[ids] == 4) & (head.elm.tag1[ids] == 2))
    assert not np.intersect1d(cortex_ids, target_ids, assume_unique=True).size
    np.testing.assert_allclose(head.nodes.node_coord[head.elm.node_number_list[target_ids] - 1].mean(axis=1),
                               target_centers, rtol=0, atol=1e-5)
    brain_ids = np.flatnonzero((head.elm.elm_type == 4) & np.isin(head.elm.tag1, [1, 2]))
    brain_vol = tetra_volumes(head, brain_ids)
    cortex_pos = np.searchsorted(brain_ids, cortex_ids)
    target_pos = np.searchsorted(brain_ids, target_ids)
    cortex_vol, target_vol = brain_vol[cortex_pos], brain_vol[target_pos]
    np.testing.assert_allclose(cortex_vol*1000, cortex_weights, rtol=2e-10, atol=1e-12)
    np.testing.assert_allclose(target_vol*1000, target_weights, rtol=2e-10, atol=1e-12)
    totals = {"cortex": cortex_vol.sum(), "brain": brain_vol.sum(), "off_target": brain_vol.sum()-target_vol.sum()}
    np.testing.assert_allclose([totals["cortex"]*1000, target_vol.sum()*1000],
                               [target["cortical_volume_mm3"], target["roi_volume_mm3"]], rtol=1e-10)
    del head
    summary, main, finely_sampled = [], [], []
    audit = {"mesh_sha256": target["mesh_sha256"], "roi_sha256": solve["roi_file_sha256"],
             "placement_sha256": solve["placement_sha256"], "native_coordinates_sha256": node_digest,
             "topology_sha256": target["topology_sha256"], "elements_cortex": len(cortex_ids),
             "elements_brain": len(brain_ids), "elements_target": len(target_ids), "fields": {}}
    for coil in prior.COILS:
        print(f"Measuring {subject} {coil}", flush=True)
        path = CACHE / subject / "fields" / f"{coil}.msh"
        mesh = simnibs.read_msh(str(path))
        assert topology_hash(mesh) == target["topology_sha256"] and coordinate_hash(mesh) == node_digest
        strength = np.linalg.norm(mesh.field["E"].value[brain_ids], axis=1)
        del mesh
        if np.any(~np.isfinite(strength)):
            raise ValueError(f"Nonfinite brain field in {path}")
        target_e, cortex_e = strength[target_pos], strength[cortex_pos]
        target_mean = np.average(target_e, weights=target_weights)
        p99 = weighted_quantile(cortex_e, cortex_weights, .99)
        old = checks[(subject, coil)]
        np.testing.assert_allclose([p99, target_mean], [old["fixed_cortex_P99_V_m"], old["fixed_target_V_m"]],
                                   rtol=1e-10, atol=1e-10)
        curves = {"cortex": Exceedance(cortex_e, cortex_vol),
                  "brain": Exceedance(strength, brain_vol), "target": Exceedance(target_e, target_vol)}
        half_volume = curves["cortex"].above(.5*p99)
        audit["fields"][coil] = {"file_sha256": sha256(path), "fixed_cortex_P99_V_m": float(p99),
                                 "fixed_target_mean_V_m": float(target_mean)}
        sources = {"cortex_P99_30": p99, "target_mean_1": target_mean}
        for condition, level in CONDITIONS.items():
            reference = sources[condition]
            gain = level / reference
            final_p99 = weighted_quantile(cortex_e * gain, cortex_weights, .99)
            final_target = np.average(target_e * gain, weights=target_weights)
            np.testing.assert_allclose([final_p99, final_target], [gain*p99, gain*target_mean], rtol=1e-13)
            if condition == "cortex_P99_30":
                np.testing.assert_allclose([final_p99, final_target],
                                           [old["matched_cortex_P99_V_m"], old["matched_target_V_m"]], rtol=1e-10)
            else:
                np.testing.assert_allclose(final_target, level, rtol=1e-13)
            summary.append(dict(subject=subject, coil=coil, condition=condition, gain=gain,
                didt_A_per_us=gain*prior.DRIVE/1e6, cortex_P99_V_m=final_p99, target_mean_V_m=final_target,
                cortex_total_cm3=totals["cortex"], brain_total_cm3=totals["brain"],
                off_target_total_cm3=totals["off_target"], cortex_half_P99_threshold_V_m=.5*final_p99,
                cortex_half_P99_cm3=half_volume, cortex_half_P99_percent=100*half_volume/totals["cortex"]))
            main.extend(volume_rows(subject, coil, condition, THRESHOLDS, reference, level, curves, totals))
            upper = int(5*np.floor(gain*strength.max()/5) + 5)
            grid = np.r_[np.arange(0, 60.5, .5), np.arange(65, upper+1, 5)]
            detail = volume_rows(subject, coil, condition, grid, reference, level, curves, totals)
            for mask in ("cortex", "brain", "off_target"):
                np.testing.assert_allclose(detail[0][mask+"_percent"], 100, atol=1e-8, rtol=0)
                assert detail[-1][mask+"_cm3"] == 0
            finely_sampled.extend(detail)
    return summary, main, finely_sampled, audit


def main():
    manifest = json.loads((COHORT / "head_sources.json").read_text())["subjects"]
    table = prior.load_csv(COHORT / "results" / "per_subject.csv")
    checks = {(r["subject"], r["coil"]): {k: float(r[k]) for k in
              ("fixed_target_V_m", "fixed_cortex_P99_V_m", "matched_cortex_P99_V_m", "matched_target_V_m")}
              for r in table if r["group"] == "replication"}
    expected = {(s, c) for s in prior.COHORT for c in prior.COILS}
    assert set(checks) == expected
    conditions, thresholds, curves, subjects = [], [], [], {}
    for subject in prior.COHORT:
        c, t, full, subjects[subject] = measure_subject(subject, checks, manifest[subject])
        conditions.extend(c)
        thresholds.extend(t)
        curves.extend(full)
    assert len(conditions) == 30 and len(thresholds) == len(conditions)*len(THRESHOLDS)
    metadata = {
        "source": "Previously solved 100 A/us quasi-static SimNIBS total element E vector fields (V/m)",
        "simnibs_version": simnibs.__version__, "source_drive_A_per_us": prior.DRIVE/1e6,
        "condition_reference_V_m": CONDITIONS, "primary_absolute_threshold_V_m": 15.0,
        "sensitivity_absolute_thresholds_V_m": THRESHOLDS.tolist(), "threshold_comparison": "inclusive >=",
        "volume_units": "cm3", "percent_definition": "100 * above-threshold volume / full mask volume",
        "masks": {"cortex": "GM tetrahedra with subject native FreeSurfer labels 3 or 42",
                  "brain": "all GM and WM tetrahedra (SimNIBS tissue tags 2 and 1)",
                  "target": "native left thalamus label 10 and GM tetrahedral barycentre within 5 mm of protocol centre",
                  "off_target": "brain mask minus target tetrahedra"},
        "curve_grid_V_m": "0..60 by 0.5; then 65,70,... to just above the condition-specific maximum brain |E|",
        "percentile": "volume-weighted empirical inverse CDF, smallest field at >=99% cumulative cortical volume",
        "shape_metric": "cortical volume at >=0.5 times its own volume-weighted cortical P99",
        "subjects": subjects,
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=RESULTS) as scratch:
        work = Path(scratch)
        for file, rows in [("conditions.csv", conditions), ("thresholds.csv", thresholds), ("curves.csv", curves)]:
            prior.save_csv(work/file, rows)
        (work/"provenance.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
        for path in work.iterdir():
            os.replace(path, RESULTS/path.name)
    print(f"Wrote {len(conditions)} conditions, {len(thresholds)} threshold rows and {len(curves)} curve rows.")


if __name__ == "__main__":
    main()
