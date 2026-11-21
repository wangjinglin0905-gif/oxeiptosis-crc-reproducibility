"""Verify patient-independent technical split counts from both frozen spatial series."""
from pathlib import Path
import numpy as np
import pandas as pd

tables = Path(__file__).resolve().parents[1] / 'tables'
reported = pd.read_csv(tables/'S43_OX4_half_direction.csv')
for cohort, section_file, split_file in [
    ('Valdeolivas', 'S33_section_results.csv', 'S37_technical_splits.csv'),
    ('GSE294385', 'S39_section_effects.csv', 'S43_split_results.csv'),
]:
    sections = pd.read_csv(tables/section_file).query("construct == 'OX4'").set_index('sample')
    splits = pd.read_csv(tables/split_file).query("construct == 'OX4'").copy()
    full = np.sign(splits['sample'].map(sections.cliff_delta))
    assert full.notna().all() and not (full == 0).any()
    splits['both'] = (np.sign(splits.delta_half_a) == full) & (np.sign(splits.delta_half_b) == full)
    aggregated = splits.groupby('sample').both.agg(['sum', 'count']).sort_index()
    expected = reported[reported.cohort == cohort].set_index('sample').sort_index()
    assert aggregated.index.equals(expected.index)
    np.testing.assert_array_equal(aggregated['sum'], expected.both_halves_agree)
    np.testing.assert_array_equal(aggregated['count'], expected.n_splits)
    print(cohort, int(aggregated['sum'].sum()), '/', int(aggregated['count'].sum()))
