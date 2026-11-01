"""
Pull NOAA GEFS analysis (observed) and GEFS 35-day forecasts from dynamical.org
and aggregate them to population-weighted state / HHS region / national daily
and weekly values.

Data:
- https://data.dynamical.org/noaa/gefs/analysis/latest.zarr          (2000-present, 3-hourly)
- https://data.dynamical.org/noaa/gefs/forecast-35-day/latest.zarr   (00Z inits 2020-10-01-present,
                                                                      31 members, 35 days)

Population weights come from the 2020 Census county population-weighted centroids:
each county's population is assigned to the nearest 0.25 degree grid cell.

Variables: 2m temperature (C), precipitable water (kg/m^2) and, from 2020 on (when
the analysis has 2m RH), absolute humidity (g/m^3) computed per grid cell.

Usage:
    python gefs_data.py analysis --start 2000-01-01 --end 2025-07-01
    python gefs_data.py forecasts --start 2020-10-01 --end 2025-06-30
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd
import xarray as xr

DATA_DIR = Path(__file__).parent / "data"
ANALYSIS_URL = "https://data.dynamical.org/noaa/gefs/analysis/latest.zarr"
FORECAST_URL = "https://data.dynamical.org/noaa/gefs/forecast-35-day/latest.zarr"
GRID_RES = 0.25

# HHS regions by state FIPS, used for ILINet region-level rows
HHS_REGIONS = {
    1: ["09", "23", "25", "33", "44", "50"],
    2: ["34", "36", "72"],
    3: ["10", "11", "24", "42", "51", "54"],
    4: ["01", "12", "13", "21", "28", "37", "45", "47"],
    5: ["17", "18", "26", "27", "39", "55"],
    6: ["05", "22", "35", "40", "48"],
    7: ["19", "20", "29", "31"],
    8: ["08", "30", "38", "46", "49", "56"],
    9: ["04", "06", "15", "32"],
    10: ["02", "16", "41", "53"],
}

ZARR_VARS = {"temp": "temperature_2m", "pw": "precipitable_water_atmosphere", "rh": "relative_humidity_2m"}
# GEFS analysis has no 2m RH before 2020, so AH is only available 2020+;
# precipitable water (kg/m^2) is the humidity variable used for features.
ANALYSIS_VARS = ("temp", "pw", "ah")
FORECAST_VARS = ("temp", "pw")

# separate bounding boxes so we don't pull the whole Pacific
REGION_OF_STATE = {"02": "AK", "15": "HI", "72": "PR"}


def load_cell_weights() -> pd.DataFrame:
    """Population per (state, grid cell) from Census 2020 county pop-weighted centroids."""
    cos = []
    for f in ["CenPop2020_Mean_CO.txt", "CenPop2020_Mean_CO72.txt"]:
        d = pd.read_csv(DATA_DIR / f, dtype={"STATEFP": str}, encoding="utf-8-sig")
        cos.append(d)
    co = pd.concat(cos, ignore_index=True)
    co["lat_i"] = np.round((90 - co["LATITUDE"]) / GRID_RES).astype(int)
    co["lon_i"] = np.round((co["LONGITUDE"] + 180) / GRID_RES).astype(int)
    cells = co.groupby(["STATEFP", "lat_i", "lon_i"], as_index=False)["POPULATION"].sum()
    cells = cells.rename(columns={"STATEFP": "location", "POPULATION": "pop"})
    cells["box"] = cells["location"].map(REGION_OF_STATE).fillna("CONUS")
    return cells


def state_utc_offsets(cells: pd.DataFrame) -> pd.Series:
    """Approximate local-time offset (hours) per state from pop-weighted mean longitude."""
    lon = cells["lon_i"] * GRID_RES - 180
    mean_lon = (lon * cells["pop"]).groupby(cells["location"]).sum() / cells.groupby("location")["pop"].sum()
    return np.round(mean_lon / 15).astype(int)


def absolute_humidity(temp_c, rh_pct):
    """Absolute humidity (g/m^3) from temperature (C) and relative humidity (%)."""
    es = 6.112 * np.exp(17.67 * temp_c / (temp_c + 243.5))
    return es * rh_pct * 2.1674 / (273.15 + temp_c)


class BoxAggregator:
    """Pull a lat/lon box from a GEFS zarr and reduce grid cells to state population-weighted means."""

    def __init__(self, cells: pd.DataFrame):
        self.cells = cells
        self.boxes = {}
        for box, c in cells.groupby("box"):
            lat0, lat1 = c["lat_i"].min(), c["lat_i"].max()
            lon0, lon1 = c["lon_i"].min(), c["lon_i"].max()
            states = sorted(c["location"].unique())
            w = np.zeros((len(c), len(states)), dtype=np.float32)
            col = c["location"].map({s: i for i, s in enumerate(states)}).to_numpy()
            w[np.arange(len(c)), col] = c["pop"].to_numpy()
            w /= w.sum(axis=0, keepdims=True)
            self.boxes[box] = dict(
                lat=slice(lat0, lat1 + 1), lon=slice(lon0, lon1 + 1),
                rows=(c["lat_i"] - lat0).to_numpy(), cols=(c["lon_i"] - lon0).to_numpy(),
                weights=w, states=states,
            )

    def reduce(self, ds: xr.Dataset, box: str, variables=("temp", "pw"), **isel) -> xr.Dataset:
        """
        Load the requested variables for one box (computing AH per cell from
        temperature and RH if asked), and return a Dataset with a 'location'
        dimension in place of latitude/longitude.
        """
        b = self.boxes[box]
        src = sorted({ZARR_VARS[v] for v in variables if v != "ah"}
                     | ({"temperature_2m", "relative_humidity_2m"} if "ah" in variables else set()))
        sub = ds[src].isel(latitude=b["lat"], longitude=b["lon"], **isel).load()
        cell = {s: sub[s].values[..., b["rows"], b["cols"]] for s in src}
        arrs = {v: cell[ZARR_VARS[v]] for v in variables if v != "ah"}
        if "ah" in variables:
            arrs["ah"] = absolute_humidity(cell["temperature_2m"], cell["relative_humidity_2m"])
        lead_dims = [d for d in sub[src[0]].dims if d not in ("latitude", "longitude")]
        coords = {d: sub[d].values for d in lead_dims}
        coords["location"] = b["states"]
        out = {}
        for name, arr in arrs.items():
            # NaN-safe weighted mean over cells (renormalize over non-missing cells)
            ok = ~np.isnan(arr)
            num = np.nan_to_num(arr) @ b["weights"]
            den = ok.astype(np.float32) @ b["weights"]
            out[name] = (lead_dims + ["location"], np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan))
        return xr.Dataset(out, coords=coords)


def add_aggregate_locations(df: pd.DataFrame, value_cols: list, by: list, cells: pd.DataFrame) -> pd.DataFrame:
    """Add HHS region and US rows as population-weighted means of state rows."""
    pop = cells.groupby("location")["pop"].sum()
    groups = {f"Region {r}": s for r, s in HHS_REGIONS.items()}
    groups["US"] = [s for s in pop.index if s != "72"]
    extra = []
    for name, states in groups.items():
        d = df[df["location"].isin(states)].copy()
        w = d["location"].map(pop)
        num = (d[value_cols].mul(w, axis=0)).groupby([d[b] for b in by]).sum()
        den = (d[value_cols].notna().mul(w, axis=0)).groupby([d[b] for b in by]).sum()
        extra.append((num / den).reset_index().assign(location=name))
    return pd.concat([df] + extra, ignore_index=True)


def to_local_daily(times: pd.DatetimeIndex, values: np.ndarray, states: list, offsets: pd.Series) -> pd.DataFrame:
    """values: (time, location) -> long daily means in local time."""
    frames = []
    for j, s in enumerate(states):
        local = times + pd.to_timedelta(int(offsets[s]), unit="h")
        frames.append(pd.DataFrame({"day": local.floor("D"), "location": s, "v": values[:, j]}))
    d = pd.concat(frames)
    return d.groupby(["location", "day"], as_index=False)["v"].mean()


def week_end(days: pd.Series) -> pd.Series:
    """Saturday ending the MMWR (Sun-Sat) week containing each day."""
    return days + pd.to_timedelta((5 - days.dt.weekday) % 7, unit="D")


@click.group()
def cli():
    pass


@cli.command()
@click.option("--start", default="2000-01-01")
@click.option("--end", default="2025-07-01")
def analysis(start, end):
    """Daily population-weighted temp and AH from GEFS analysis -> data/gefs_analysis_daily.csv"""
    cells = load_cell_weights()
    offsets = state_utc_offsets(cells)
    agg = BoxAggregator(cells)
    ds = xr.open_zarr(ANALYSIS_URL, decode_timedelta=True)
    out_path = DATA_DIR / "gefs_analysis_daily.csv"
    done = set()
    if out_path.exists():
        done = set(pd.read_csv(out_path, usecols=["block"])["block"].unique())
    i0, i1 = np.searchsorted(ds.time.values, [np.datetime64(start), np.datetime64(end)])
    block = 1440  # matches zarr time chunking
    for c in range(i0 // block * block, i1, block):
        sl = slice(max(c, i0), min(c + block, i1))
        tag = str(pd.Timestamp(ds.time.values[sl.start]).date())
        if tag in done:
            continue
        frames = []
        for box in agg.boxes:
            red = agg.reduce(ds, box, variables=ANALYSIS_VARS, time=sl)
            times = pd.DatetimeIndex(red["time"].values)
            states = list(red["location"].values)
            d = [to_local_daily(times, red[v].values, states, offsets).set_index(["location", "day"])["v"].rename(v)
                 for v in ANALYSIS_VARS]
            frames.append(pd.concat(d, axis=1))
        daily = pd.concat(frames).reset_index()
        daily["block"] = tag
        daily.to_csv(out_path, mode="a", header=not out_path.exists(), index=False)
        print(f"analysis block {tag}: {len(daily)} rows", flush=True)


@cli.command()
@click.option("--start", default="2020-10-01")
@click.option("--end", default="2025-06-30")
@click.option("--months", default="9,10,11,12,1,2,3,4,5,6", help="init months to pull")
def forecasts(start, end, months):
    """
    Weekly ensemble mean/sd of population-weighted temp and AH for MMWR weeks
    1-3 after each Wednesday 00Z init -> data/gefs_forecast_weekly.csv

    Week k ends on (init Wednesday + 3 days) + 7*(k-1) days, i.e. k=1 is the week
    containing the Wednesday (only Wed-Sat forecast days are kept; Sun-Tue are
    filled from the analysis later), k=2 and k=3 are fully forecast.
    """
    cells = load_cell_weights()
    offsets = state_utc_offsets(cells)
    agg = BoxAggregator(cells)
    ds = xr.open_zarr(FORECAST_URL, decode_timedelta=True)
    months = {int(m) for m in months.split(",")}
    inits = [d for d in pd.date_range(start, end, freq="W-WED") if d.month in months]
    out_path = DATA_DIR / "gefs_forecast_weekly.csv"
    done = set()
    if out_path.exists():
        done = set(pd.read_csv(out_path, usecols=["init_date"])["init_date"].unique())
    lead_hours = ds.lead_time.values / np.timedelta64(1, "h")
    n_lead = int(np.searchsorted(lead_hours, 19 * 24)) + 1
    for init in inits:
        tag = str(init.date())
        if tag in done:
            continue
        if np.datetime64(init) not in ds.init_time.values:
            print(f"missing init {tag}", flush=True)
            continue
        i_init = int(np.where(ds.init_time.values == np.datetime64(init))[0][0])
        sat = init + pd.Timedelta(days=3)
        frames = []
        for box in agg.boxes:
            red = agg.reduce(ds, box, variables=FORECAST_VARS, init_time=i_init, lead_time=slice(0, n_lead))
            times = init + pd.to_timedelta(red["lead_time"].values)
            for var in FORECAST_VARS:
                arr = red[var].values  # (member, lead, location)
                for j, s in enumerate(red["location"].values):
                    local = times + pd.to_timedelta(int(offsets[s]), unit="h")
                    d = pd.DataFrame(arr[:, :, j].T, index=local.floor("D"))
                    daily = d.groupby(level=0).mean()
                    daily = daily[daily.index >= init]  # local Wednesday onward
                    wk = week_end(pd.Series(daily.index, index=daily.index))
                    weekly = daily.groupby(wk.values).mean()
                    ndays = daily.groupby(wk.values).size()
                    for k in (1, 2, 3):
                        we = sat + pd.Timedelta(weeks=k - 1)
                        if we not in weekly.index:
                            continue
                        need = 4 if k == 1 else 7
                        if ndays[we] < need:
                            continue
                        m = weekly.loc[we].to_numpy()
                        frames.append(dict(init_date=tag, location=s, wk_end_date=str(we.date()), k=k,
                                           var=var, mean=float(np.mean(m)), sd=float(np.std(m, ddof=1))))
        out = pd.DataFrame(frames)
        out.to_csv(out_path, mode="a", header=not out_path.exists(), index=False)
        print(f"forecast init {tag}: {len(out)} rows", flush=True)


if __name__ == "__main__":
    cli()
