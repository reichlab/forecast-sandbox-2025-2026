"""
Real-time ("as published") state-level ILINet x percent positive for the development seasons 2018/19 and 2019/20,
from the Delphi Epidata API (https://cmu-delphi.github.io/delphi-epidata/).

iddata's ILINet series is the unweighted ILI of each state (weighted ILI for the nation) times the WHO/NREVSS
clinical-lab percent positive, divided by 100, using the current (final) values. Delphi keeps the values as first
published and at every later weekly release ("issue"). State ILI vintages start with issue 201740 and clinical-lab
vintages with issue 201839, which covers both development seasons.

Fetching. For each lag L = 0, ..., MAX_LAG - 1 and each of the two endpoints (`fluview`, `fluview_clinical`), one API
call returns, for every region and epiweek, the value published L weeks after the epiweek ended (issue = epiweek + L).
That is 2 * MAX_LAG calls in total, well under the anonymous limit of 60 per hour. The results are cached in
analysis/peak-models/ilinet-validation/epidata/.

As-of reconstruction. At a reference date r (a Saturday), the latest published epiweek is the one ending r - 7, whose
issue number we call X. The value of epiweek w as known at r is its value at lag X - w if X - w < MAX_LAG; weeks with
X - w >= MAX_LAG are treated as final and take iddata's value. (Revisions beyond lag 4 are under 1% at the median.)
If a lag is missing for some week (occasional late or skipped releases), the nearest earlier available lag is used; if
none is available, the week is missing at r.

Regions: the 50 states other than Florida (no ILINet), DC, and New York as `ny_minus_jfk` (NY excluding NYC, which
matches iddata's "New York": NYC has no clinical-lab positivity, so iddata's average reduces to upstate NY), plus the
nation (`nat`, using weighted ILI). Puerto Rico and the Virgin Islands are excluded (see validate_ilinet.py).

Usage (from the repository root):
    python analysis/peak-models/epidata_ilinet.py            # fetch (resumable) and summarize
"""
import datetime
import json
import sys
import time
from pathlib import Path

import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
CACHE = HERE / "ilinet-validation" / "epidata"
BASE = "https://api.delphi.cmu.edu/epidata/"
MAX_LAG = 10
EPIWEEKS = "201740-202035"

STATE_FIPS = {
    "al": "01", "ak": "02", "az": "04", "ar": "05", "ca": "06", "co": "08", "ct": "09", "de": "10", "dc": "11",
    "ga": "13", "hi": "15", "id": "16", "il": "17", "in": "18", "ia": "19", "ks": "20", "ky": "21", "la": "22",
    "me": "23", "md": "24", "ma": "25", "mi": "26", "mn": "27", "ms": "28", "mo": "29", "mt": "30", "ne": "31",
    "nv": "32", "nh": "33", "nj": "34", "nm": "35", "ny_minus_jfk": "36", "nc": "37", "nd": "38", "oh": "39",
    "ok": "40", "or": "41", "pa": "42", "ri": "44", "sc": "45", "sd": "46", "tn": "47", "tx": "48", "ut": "49",
    "vt": "50", "va": "51", "wa": "53", "wv": "54", "wi": "55", "wy": "56", "nat": "US",
}


def _get(endpoint: str, lag: int | None, epiweeks: str = EPIWEEKS, tag: str = "") -> list[dict]:
    """One Epidata call for all regions at one lag (lag None: the latest values), cached as JSON."""
    path = CACHE / f"{endpoint}{tag}_lag{lag if lag is not None else 'latest'}.json"
    if path.exists():
        return json.loads(path.read_text())
    params = {"regions": ",".join(STATE_FIPS), "epiweeks": epiweeks}
    if lag is not None:
        params["lag"] = lag
    url = BASE + endpoint + "/?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                remaining, j = resp.headers.get("x-my-remaining"), json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code != 429:
                raise
            wait = int(e.headers.get("x-my-reset", 600)) + 5
            print(f"rate limited; waiting {wait}s", flush=True)
            time.sleep(wait)
            continue
        if j.get("result") not in (1, -2):  # 1: ok, -2: no results; 2 would mean truncated
            raise RuntimeError(f"{endpoint} lag {lag}: result {j.get('result')} {j.get('message')}")
        rows = j.get("epidata") or []
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(rows))
        print(f"{endpoint} lag {lag}: {len(rows)} rows (remaining calls this hour: {remaining})", flush=True)
        return rows
    raise RuntimeError(f"{endpoint} lag {lag}: gave up after repeated rate limiting")


def fetch() -> None:
    for endpoint in ("fluview", "fluview_clinical"):
        for lag in range(MAX_LAG):
            _get(endpoint, lag)


def epiweek_to_saturday(ew: int) -> datetime.date:
    """The Saturday ending MMWR epiweek `ew` (YYYYWW)."""
    year, week = divmod(int(ew), 100)
    jan4 = datetime.date(year, 1, 4)
    # MMWR week 1 is the Sunday-Saturday week containing January 4
    first_sunday = jan4 - datetime.timedelta(days=(jan4.weekday() + 1) % 7)
    return first_sunday + datetime.timedelta(weeks=week - 1, days=6)


def lag_table() -> pd.DataFrame:
    """Long table: location (FIPS), wk_end_date, lag, ili (unweighted for states, weighted for US), percent_positive,
    and clinical-lab influenza A and B positives (total_a, total_b)."""
    frames = []
    for lag in range(MAX_LAG):
        ili = pd.DataFrame(_get("fluview", lag))
        if len(ili):
            ili["ili_value"] = np.where(ili["region"] == "nat", ili["wili"], ili["ili"])
            ili = ili[["region", "epiweek", "ili_value"]]
        clin = pd.DataFrame(_get("fluview_clinical", lag))
        if len(clin):
            clin = clin[["region", "epiweek", "percent_positive", "total_a", "total_b"]]
        df = ili.merge(clin, on=["region", "epiweek"], how="outer").assign(lag=lag)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["location"] = out["region"].map(STATE_FIPS)
    out["wk_end_date"] = pd.to_datetime([epiweek_to_saturday(e) for e in out["epiweek"]])
    return out.drop(columns="region")


def as_of_values(table: pd.DataFrame, ref_date: datetime.date, final: pd.DataFrame) -> pd.DataFrame:
    """
    ILINet x percent positive (inc) as known at ref_date, per location and wk_end_date, for weeks ending by ref_date - 7.
    `final` holds iddata's final values (location, wk_end_date, inc) and is used for weeks at least MAX_LAG weeks old.
    ILI and positivity are each taken at the largest available lag not exceeding the week's age at ref_date.
    """
    last = pd.Timestamp(ref_date - datetime.timedelta(days=7))
    t = table.loc[table["wk_end_date"] <= last].copy()
    t["age"] = ((last - t["wk_end_date"]).dt.days // 7).astype(int)
    recent = t.loc[t["age"] < MAX_LAG]
    parts = {}
    for col in ("ili_value", "percent_positive"):
        v = recent.loc[recent[col].notna() & (recent["lag"] <= recent["age"])]
        v = v.sort_values("lag").groupby(["location", "wk_end_date"])[col].last()  # largest lag <= age
        parts[col] = v
    rt = pd.concat(parts, axis=1).reset_index()
    rt["inc"] = rt["ili_value"] * rt["percent_positive"] / 100.0
    rt = rt[["location", "wk_end_date", "inc"]]
    old = final.loc[(final["wk_end_date"] <= last - pd.Timedelta(weeks=MAX_LAG)), ["location", "wk_end_date", "inc"]]
    return pd.concat([old, rt], ignore_index=True).sort_values(["location", "wk_end_date"]).reset_index(drop=True)


NHSN_EPIWEEKS = "202235-202535"


def nhsn_clinical_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Clinical-lab A and B positives for the NHSN test seasons (epiweeks NHSN_EPIWEEKS): a lag table (lags 0..MAX_LAG-1,
    as lag_table) and the latest values. One API call per lag plus one for the latest values.
    """
    frames = []
    for lag in list(range(MAX_LAG)) + [None]:
        c = pd.DataFrame(_get("fluview_clinical", lag, epiweeks=NHSN_EPIWEEKS, tag="_nhsn"))
        c = c[["region", "epiweek", "total_a", "total_b"]].assign(lag=lag if lag is not None else -1)
        frames.append(c)
    out = pd.concat(frames, ignore_index=True)
    out["location"] = out["region"].map(STATE_FIPS)
    out["wk_end_date"] = pd.to_datetime([epiweek_to_saturday(e) for e in out["epiweek"]])
    out = out.drop(columns="region")
    return out.loc[out["lag"] >= 0].reset_index(drop=True), out.loc[out["lag"] < 0].drop(columns="lag")


def strain_as_of(table: pd.DataFrame, ref_date: datetime.date, season: str, week_map: pd.DataFrame,
                 final: dict) -> dict:
    """
    Type/subtype counts for `season` as known at ref_date, in the format of idmodels.peak.extra_features.type_arrays
    ({(geography, season): array (4, 53) of A, B, A(H1), A(H3) by season week}), starting from `final` (final counts
    for all seasons). For states and the nation, A and B for weeks published by ref_date are replaced by their
    as-published clinical-lab values (largest lag available, weeks at least MAX_LAG old keep final values) and
    later weeks are removed; A(H1) and A(H3) keep their final values (no subtype vintages are available), which
    uses revisions to subtype counts that were not yet known.
    """
    last = pd.Timestamp(ref_date - datetime.timedelta(days=7))
    t = table.loc[(table["wk_end_date"] <= last) & table["total_a"].notna()].copy()
    t["age"] = ((last - t["wk_end_date"]).dt.days // 7).astype(int)
    t = t.loc[(t["age"] < MAX_LAG) & (t["lag"] <= t["age"])]
    t = t.sort_values("lag").groupby(["location", "wk_end_date"])[["total_a", "total_b"]].last().reset_index()
    t = t.merge(week_map, on="wk_end_date")
    out = dict(final)
    last_sw = int(week_map.loc[week_map["wk_end_date"] == last, "season_week"].iloc[0])
    for geo in {g for g, s in final if s == season} | set(t["location"].dropna()):
        arr = final.get((geo, season), np.full((4, 53), np.nan)).copy()
        arr[:, last_sw:] = np.nan  # nothing after the last published week
        rows = t.loc[t["location"] == geo]
        if len(rows):
            wk = rows["season_week"].to_numpy().astype(int) - 1
            arr[0, wk], arr[1, wk] = rows["total_a"].to_numpy(float), rows["total_b"].to_numpy(float)
        out[(geo, season)] = arr
    return out


def revision_vintages(table: pd.DataFrame, ref_date: datetime.date) -> pd.DataFrame:
    """
    ILINet x positivity vintages known at ref_date, in the format RevisionModel.fit expects (location, wk_end_date,
    inc, as_of). Each release X (identified by its latest epiweek) contributes the values it published for weeks at
    lags 0..MAX_LAG-1 (both ILI and positivity at that lag), for locations with a lag-0 value in that release. The
    "final" vintage (as_of = ref_date) holds each week's value at the largest lag available by ref_date (at most
    MAX_LAG - 1), so revision ratios are measured against what was known at ref_date and nothing later is used.
    """
    last = pd.Timestamp(ref_date - datetime.timedelta(days=7))
    t = table.dropna(subset=["ili_value", "percent_positive"]).copy()
    t["as_of"] = t["wk_end_date"] + pd.to_timedelta(7 * t["lag"], unit="D")  # the release's latest week
    t = t.loc[t["as_of"] <= last]
    t["inc"] = t["ili_value"] * t["percent_positive"] / 100.0
    has_lag0 = t.loc[t["lag"] == 0, ["as_of", "location"]].drop_duplicates()
    releases = t.merge(has_lag0, on=["as_of", "location"])[["location", "wk_end_date", "inc", "as_of"]]
    final = (
        t.sort_values("lag").groupby(["location", "wk_end_date"])["inc"].last().reset_index().assign(as_of=last)
    )
    return pd.concat([releases, final], ignore_index=True)


def main():
    fetch()
    table = lag_table()
    print(table.groupby("lag")[["ili_value", "percent_positive"]].count().to_string())
    sys.exit(0)


if __name__ == "__main__":
    main()
