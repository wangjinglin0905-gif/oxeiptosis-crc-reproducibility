"""Reproduce the frozen OX4 morphology endpoint and count-split readouts."""

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from threadpoolctl import threadpool_limits

from prepare_inputs import OUT, SAMPLES, SEED, SETS, ZONES, sha, write_json


def cliff_delta(scores, tumor, stroma):
    use = tumor | stroma
    n_tumor, n_stroma = int(tumor.sum()), int(stroma.sum())
    ranked = stats.rankdata(scores[use], axis=0)
    u = ranked[tumor[use]].sum(axis=0) - n_tumor * (n_tumor + 1) / 2
    return 2 * u / (n_tumor * n_stroma) - 1


def spearman_columns(a, b):
    ra, rb = stats.rankdata(a, axis=0), stats.rankdata(b, axis=0)
    ra -= ra.mean(axis=0)
    rb -= rb.mean(axis=0)
    denominator = np.sqrt((ra**2).sum(axis=0) * (rb**2).sum(axis=0))
    return np.divide((ra * rb).sum(axis=0), denominator,
                     out=np.full(a.shape[1], np.nan), where=denominator > 0)


def compute_scores(counts, library, center, sd, weights):
    scaled = (np.log1p(10000 * counts / library[:, None]) - center) / sd
    return scaled @ weights


def subgroup_delta(scores, tumor, stroma):
    if min(tumor.sum(), stroma.sum()) < 30:
        return np.full(scores.shape[1], np.nan)
    return cliff_delta(scores, tumor, stroma)


def verify_delta(x, tumor, stroma):
    sorted_stroma = np.sort(x[stroma])
    wins = np.searchsorted(sorted_stroma, x[tumor], side="left").sum()
    losses = stroma.sum() * tumor.sum() - np.searchsorted(
        sorted_stroma, x[tumor], side="right").sum()
    return (wins - losses) / (tumor.sum() * stroma.sum())


def section(gate):
    sample = gate["sample"]
    directory = OUT / "sections" / sample
    checkpoint = directory / "eligible_counts.npz"
    with np.load(checkpoint) as record:
        counts, library, center, sd, weights, names, labels, zones, barcodes = [
            record[key] for key in ("count", "library", "center", "sd", "weights", "names",
                                    "compartments", "zones", "barcodes")
        ]
    tumor, stroma = labels == "tumor", labels == "stroma"
    assert (len(library), tumor.sum(), stroma.sum()) == (
        gate["eligible_spots"], gate["n_tumor"], gate["n_stroma"])
    fixed = np.array([label in SETS for label in names])
    rest = library - counts.sum(axis=1)
    assert (rest >= 0).all()
    log_count = np.log1p(library)
    unscaled = compute_scores(counts, library, center, sd, weights)
    composite_sd = unscaled.std(axis=0, ddof=1)
    assert (composite_sd > 0).all()
    scores = unscaled / composite_sd
    assert np.max(abs(scores.mean(axis=0))) < 1e-9
    delta = cliff_delta(scores, tumor, stroma)
    validation_delta = [verify_delta(x, tumor, stroma) for x in scores.T]
    formula_error = float(max(abs(delta - validation_delta)))
    assert formula_error < 1e-12

    design = np.column_stack([np.ones(len(log_count)), log_count])
    residual = scores - design @ np.linalg.lstsq(design, scores, rcond=None)[0]
    away = zones != ZONES[-1]
    zone_effects = []
    for zone in ZONES:
        loc = zones == zone
        estimate = subgroup_delta(scores, tumor & loc, stroma & loc)
        zone_effects.append(estimate)
    wall = np.stack(zone_effects[:3])
    wall_n = int(np.isfinite(wall[:, 0]).sum())
    wall_average = np.nansum(wall, axis=0) / wall_n if wall_n else np.full(len(names), np.nan)
    mean_difference = scores[tumor].mean(axis=0) - scores[stroma].mean(axis=0)

    repeats = []
    max_split_error = 0.0
    for i in range(30):
        rng = np.random.default_rng(np.random.SeedSequence([SEED, list(SAMPLES).index(sample), i]))
        left = rng.binomial(counts, .5)
        left_library = left.sum(axis=1) + rng.binomial(rest, .5)
        right_library = library - left_library
        assert (left_library > 0).all() and (right_library > 0).all()
        a = compute_scores(left, left_library, center, sd, weights) / composite_sd
        b = compute_scores(counts - left, right_library, center, sd, weights) / composite_sd
        rho = spearman_columns(a, b)
        rho_tumor = spearman_columns(a[tumor], b[tumor])
        rho_stroma = spearman_columns(a[stroma], b[stroma])
        da, db = cliff_delta(a, tumor, stroma), cliff_delta(b, tumor, stroma)
        if i == 0:
            for j in np.flatnonzero(fixed):
                expected = [stats.spearmanr(a[g, j], b[g, j]).statistic
                            for g in (slice(None), tumor, stroma)]
                observed = [rho[j], rho_tumor[j], rho_stroma[j]]
                assert np.array_equal(np.isnan(expected), np.isnan(observed))
                for e, o in zip(expected, observed):
                    if np.isfinite(e):
                        max_split_error = max(max_split_error, abs(e - o))
            assert max_split_error < 1e-12
        repeats.append(pd.DataFrame(dict(construct=names, sample=sample, rep=i, rho_all=rho,
                                         rho_tumor=rho_tumor, rho_stroma=rho_stroma,
                                         delta_half_a=da, delta_half_b=db,
                                         mean_absolute_delta_shift=(abs(da-delta)+abs(db-delta))/2)))
    repeat = pd.concat(repeats, ignore_index=True)
    repeat.to_csv(directory / "split_results.csv", index=False)
    med = repeat.groupby("construct").median(numeric_only=True).loc[names]
    summary = pd.DataFrame(dict(sample=sample, patient=gate["patient"], construct=names,
                                n_eligible=len(labels), n_tumor=tumor.sum(), n_stroma=stroma.sum(),
                                cliff_delta=delta, mean_difference_sd=mean_difference,
                                depth_adjusted_delta=cliff_delta(residual, tumor, stroma),
                                n_nonborder_tumor=int((tumor & away).sum()),
                                n_nonborder_stroma=int((stroma & away).sum()),
                                nonborder_delta=subgroup_delta(scores, tumor & away, stroma & away),
                                eligible_wall_strata=wall_n, wall_stratified_delta=wall_average,
                                median_split_rho=med.rho_all.to_numpy(),
                                median_tumor_rho=med.rho_tumor.to_numpy(),
                                median_stroma_rho=med.rho_stroma.to_numpy(),
                                median_half_delta_shift=med.mean_absolute_delta_shift.to_numpy(),
                                n_undefined_whole_splits=repeat["rho_all"].isna().groupby(
                                    repeat.construct).sum().loc[names].to_numpy()))
    summary.to_csv(directory / "all_construct_results.csv", index=False)
    pd.concat([pd.DataFrame(dict(sample=sample, construct=names[fixed], zone=zone,
                                n_tumor=int((tumor & (zones == zone)).sum()),
                                n_stroma=int((stroma & (zones == zone)).sum()),
                                cliff_delta=effect[fixed]))
               for zone, effect in zip(ZONES, zone_effects)], ignore_index=True).to_csv(
                   directory / "zone_effects.csv", index=False)
    scored_spots = pd.DataFrame(dict(barcode=barcodes, library=library, compartment=labels, zone=zones))
    for i in np.flatnonzero(fixed):
        scored_spots[names[i]] = scores[:, i]
    scored_spots.to_csv(directory / "fixed_scores.csv.gz", compression="gzip", index=False)
    check = dict(sample=sample, cliff_delta_error=formula_error, spearman_first_split_error=max_split_error,
                 count_checkpoint_sha256=sha(checkpoint), n_technical_rows=len(repeat))
    write_json(directory / "result_check.json", check)
    print(sample, f"delta={float(delta[0]):.4f}", f"split={med.rho_all.iloc[0]:.4f}",
          f"{len(repeat)} split rows", flush=True)
    return summary, check


def infer_patient(construct, data):
    effects = data.cliff_delta.to_numpy()
    count = len(effects)
    rng = np.random.default_rng(SEED + list(SETS).index(construct))
    boot = effects[rng.integers(0, count, size=(10000, count))].mean(axis=1)
    signs = np.array(list(itertools.product([-1, 1], repeat=count)))
    point = effects.mean()
    loo = np.array([np.delete(effects, i).mean() for i in range(count)]) if count > 1 else [np.nan]
    return dict(construct=construct, n_patients=count, mean_delta=point,
                bootstrap_lo=np.quantile(boot, .025), bootstrap_hi=np.quantile(boot, .975),
                exact_signflip_p=np.mean(abs((signs * effects).mean(axis=1)) >= abs(point) - 1e-12),
                positive_patients=int((effects > 0).sum()), negative_patients=int((effects < 0).sum()),
                leave_one_patient_out_min=min(loo), leave_one_patient_out_max=max(loo))


def summarize(sections):
    sections.to_csv(OUT / "all_construct_results.csv", index=False)
    fixed = sections[sections.construct.isin(SETS)].copy()
    fixed.to_csv(OUT / "section_effects.csv", index=False)
    cols = ["cliff_delta", "mean_difference_sd", "depth_adjusted_delta", "nonborder_delta",
            "wall_stratified_delta", "median_split_rho", "median_tumor_rho", "median_stroma_rho",
            "median_half_delta_shift"]
    patients = []
    for (patient, construct), data in fixed.groupby(["patient", "construct"]):
        effects = {name: data[name].mean(skipna=False) for name in cols}
        patients.append(dict(patient=int(patient), construct=construct, sections=len(data),
                             complete_planned_sections=len(data) == sum(
                                 p == patient for p, _ in SAMPLES.values()),
                             concordant_section_direction=(np.sign(data.cliff_delta).nunique() == 1)
                             if len(data) > 1 else None, **effects))
    patient = pd.DataFrame(patients)
    patient.to_csv(OUT / "patient_effects.csv", index=False)
    inference = pd.DataFrame([infer_patient(name, sub)
                              for name in SETS if len(sub := patient[patient.construct.eq(name)])])
    inference["BH_q"] = stats.false_discovery_control(inference.exact_signflip_p)
    inference.to_csv(OUT / "patient_inference.csv", index=False)

    background = []
    for sample, rows in sections.groupby("sample"):
        ox = rows.loc[rows.construct.eq("OX4")].iloc[0]
        matches = rows[rows.construct.str.startswith("random_")]
        if len(matches):
            assert len(matches) == 100 and matches.median_split_rho.notna().all()
            background.append(dict(sample=sample, matched_sets=len(matches),
                                   absolute_delta_percentile=100*np.mean(
                                       abs(matches.cliff_delta) <= abs(ox.cliff_delta)),
                                   split_rho_percentile=100*np.mean(
                                       matches.median_split_rho <= ox.median_split_rho),
                                   n_both_at_least_delta_and_no_greater_rho=int((
                                       abs(matches.cliff_delta).ge(abs(ox.cliff_delta)) &
                                       matches.median_split_rho.le(ox.median_split_rho)).sum())))
    if background:
        pd.DataFrame(background).to_csv(OUT / "matched_background_positions.csv", index=False)
    print(inference.to_string(index=False), flush=True)


def main():
    gates = json.loads((OUT / "input_gate.json").read_text(encoding="utf-8"))
    assert [g["sample"] for g in gates] == list(SAMPLES)
    estimates, checks = [], []
    with threadpool_limits(limits=4):
        for gate in gates:
            if gate["comparable"]:
                result, check = section(gate)
                estimates.append(result)
                checks.append(check)
    summarize(pd.concat(estimates, ignore_index=True))
    files = ["analysis_contract.md", "input_gate.json", "input_provenance.json",
             "patient_effects.csv", "patient_inference.csv", "section_effects.csv"]
    write_json(OUT / "analysis_checks.json", dict(
        completed_sections=len(checks), planned_sections=len(gates), sections=checks,
        software=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        script_sha256={name: sha(Path(__file__).with_name(name)) for name in
                       ("prepare_inputs.py", "export_annotations.R", "replicate_endpoint.py")},
        artifact_sha256={name: sha(OUT / name) for name in files}))


if __name__ == "__main__":
    main()
