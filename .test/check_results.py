"""Check that the workflow recovers the differential expression planted by
make_test_data.py (genes 1-8 up, 9-16 down in TREAT)."""
import sys

import pandas as pd

res = pd.read_csv(".test/results/diffexp/TREAT_vs_CTRL.results.tsv", sep="\t")
truth_up = {f"GENE{i:03d}" for i in range(1, 9)}
truth_down = {f"GENE{i:03d}" for i in range(9, 17)}
up = set(res.gene_id[res.status == "Up"])
down = set(res.gene_id[res.status == "Down"])

recall = (len(up & truth_up) + len(down & truth_down)) / 16
false_calls = len(up - truth_up) + len(down - truth_down)
print(f"recall = {recall:.2f}, false calls = {false_calls}")
if recall < 0.9 or false_calls > 1:
    sys.exit("Differential expression results do not match the simulated truth")
