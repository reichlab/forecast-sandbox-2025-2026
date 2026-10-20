"""
Save the Delphi Epidata archives of NHSN flu admissions and hospitals reporting (every published value, with its
publish date) to OUT_DIR, for main.py --archive_dir. Saving them once keeps a batch of runs from fetching them in every
job, and makes the runs reproducible.

  python fetch_archives.py archives
"""

import sys
from pathlib import Path

from idmodels.nowcast.delphi import nhsn_admissions_archive, nhsn_hosp_reporting_archive

out_dir = Path(sys.argv[1])
out_dir.mkdir(parents=True, exist_ok=True)
for name, archive in [("nhsn_admissions", nhsn_admissions_archive()), ("nhsn_hosp_reporting", nhsn_hosp_reporting_archive())]:
    archive.to_parquet(out_dir / f"{name}.parquet", index=False)
    print(f"{name}: {len(archive)} rows, published {archive['report_date'].min().date()} to {archive['report_date'].max().date()}")
