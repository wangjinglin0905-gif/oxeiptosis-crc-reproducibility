"""Check the GSE294385 summary from the packaged section and split tables."""

from pathlib import Path
import itertools

import numpy as np
import pandas as pd


def main():
    directory = Path(__file__).resolve().parents[1] / "tables"
    sections = pd.read_csv(directory / "S39_section_effects.csv")
    patients = pd.read_csv(directory / "S38_patient_effects.csv").set_index(["patient","construct"])
    summary = pd.read_csv(directory / "S38_patient_inference.csv").set_index("construct")
    splits = pd.read_csv(directory / "S43_split_results.csv")
    assert len(splits) == 24480
    expected = sections.set_index(["sample","construct"])["median_split_rho"].sort_index()
    observed = splits.groupby(["sample","construct"]).rho_all.median().loc[expected.index]
    np.testing.assert_allclose(expected, observed, atol=1e-12, rtol=0)
    vals = sections.groupby(["patient","construct"]).cliff_delta.mean()
    np.testing.assert_allclose(vals, patients.cliff_delta.loc[vals.index], atol=1e-12, rtol=0)
    for name, offset in [("OX4",0),("stromal",2)]:
        x = vals.xs(name, level="construct").to_numpy()
        assert len(x) == 5
        np.testing.assert_allclose(np.mean(x), summary.loc[name,"mean_delta"], atol=1e-12, rtol=0)
        rng = np.random.default_rng(20261118+offset)
        interval = np.quantile(x[rng.integers(0,5,size=(10000,5))].mean(axis=1), [.025,.975])
        np.testing.assert_allclose(interval, summary.loc[name,["bootstrap_lo","bootstrap_hi"]].astype(float),
                                   atol=1e-12,rtol=0)
        signs = np.array(list(itertools.product([-1,1],repeat=5)))
        p = np.mean(np.abs((signs*x).mean(axis=1)) >= abs(x.mean())-1e-12)
        assert p == summary.loc[name,"exact_signflip_p"]
    print("Verified 8 sections, 5 patients, 24480 split records and patient inference.")


if __name__ == "__main__":
    main()
