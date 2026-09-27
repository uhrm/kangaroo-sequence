"""Write the decade transition table as CSV.

One row per decade k in [K, K + P) and input (the entry 10^k + u or the start
c * 10^k). Decade k >= K uses the row for K + (k - K) mod P.

Run with: uv run python analysis/export_table.py [base] [output.csv]
"""

import csv
import sys

from kangaroo_sequence.decades import decade_period, transition_table

base = int(sys.argv[1]) if len(sys.argv) > 1 else 10
path = sys.argv[2] if len(sys.argv) > 2 else f"analysis/output/decade_table_base{base}.csv"
K, P = decade_period(base)

with open(path, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["k", "input", "value", "exit_u", "landmine_r", "branch_points"])
    for k, row in zip(range(K, K + P), transition_table(base)):
        for kind, inputs in (("entry", row.entries), ("start", row.starts)):
            for value, passage in inputs.items():
                branch_points = " ".join(str(d) for d, _, _ in passage.branch_points)
                writer.writerow([k, kind, value, passage.exit, passage.landmine, branch_points])
print(f"wrote {path}: decades {K}..{K + P - 1} (period {P})")
