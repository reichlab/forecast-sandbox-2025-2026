"""WISAR6_fourthroot_adaptive_t: WIS-loss AR(6) + Fourier seasonality,
pooled across locations -- fourth-root transform, learned (data-adaptive)
prior scales on phi/sigma/Fourier coefficients (inherited from
../WISAR6_fourthroot_adaptive/ unchanged), PLUS a heavier-tailed Student-t
innovation distribution in place of Gaussian. No differencing (see below
for why this is a deliberately isolated test of Student-t alone).

A from-scratch PyMC/PyTensor implementation (no idmodels/sarix dependency),
developed and validated in https://github.com/<user>/reading-literature's
pymc-tutorial/ notebooks:
  - wis_vs_likelihood_ar.py       (AR(2)/AR(6) simulation study)
  - ar2_seasonal_flusight.py      (single-location, real data)
  - ar6_fourier_pooled_flusight.py (4 locations, pooling toggle)
  - ar6_fourier_all_locations.py (all 53 locations, extended training window)

This module is a sibling of ../WISAR6_fourthroot_adaptive/model.py (which
is itself a sibling of ../WISAR6_fourthroot/, differing only in the
phi/sigma/Fourier prior structure -- see that module's docstring for the
WIS-loss/Gibbs-posterior background and the learned-prior rationale, both
of which apply unchanged here). This module adds exactly one further
change: a Student-t innovation distribution (see "Student-t innovations"
below). It is a sibling of ../WISAR6_fourthroot_adaptive_diff/, which
bundled this same Student-t change together with ordinary
first-differencing (d=1) -- that combined test found the overall WIS got
substantially worse, but the methodology could not separate whether
differencing or Student-t (or their interaction) was responsible. This
module exists specifically to isolate Student-t alone: same model as
../WISAR6_fourthroot_adaptive/ in every other respect (no differencing),
so any WIS change here is attributable to the innovation distribution
alone. See model-output/UMass-WISAR6_fourthroot/PERFORMANCE_NOTES.md's
"Attempt 2" and "Attempt 3" sections for the full diagnosis chain.

Student-t innovations: selected via `innovation_dist="studentt"` (the
default here) in fit_model/simulate_forecast_draws, in place of
"gaussian". See fit_model's and studentt_wis_loss's docstrings for the
important caveat that the degrees-of-freedom parameter `t_dof` must be
fixed, not sampled, under the WIS-loss objective (no implemented gradient
through PyTensor's betaincinv w.r.t. degrees of freedom).

This module contains the reusable model-fitting and forecasting logic;
`main.py` is the CLI entrypoint that wires it to the hub's expected
--today_date interface and hubverse output format.

Model equation, per location ell:

    z_{ell,t} = mu0_ell + S_ell(t)
                + sum_{j=1}^{6} phi_j,ell * (z_{ell,t-j} - mu0_ell - S_ell(t-j))
                + eps_{ell,t},   eps_{ell,t} ~ Normal(0, sigma_ell^2)

    S_ell(t) = sum_{k=1}^{K} [a_k,ell * sin(2*pi*k*doy(t)/365.25)
                               + b_k,ell * cos(2*pi*k*doy(t)/365.25)]

Data transform: rate-per-100k, then fourth-root, then center/scale --

    rate_t = y_t / population * 100,000
    r_t = (rate_t + offset)^0.25
    z_t = r_t / q95 - mean_scaled

where q95 and mean_scaled are the 95th percentile and mean of r_t restricted
to in-season weeks, computed once from training data, and offset matches
idmodels.constants.POWER_TRANSFORM_OFFSET (0.01). Predictions on the
transformed scale are converted back to predictions on the original scale
by inverting these transformations.

Two fitting objectives, selected via `objective`:
  - "likelihood": standard Gaussian likelihood (pm.Normal(..., observed=...))
  - "wis": Gibbs-posterior loss built from the Weighted Interval Score at the
    full 23-level FluSight quantile grid, added via pm.Potential. See
    ../../../reading-literature/scoring_estimation/generalized_bayesian_scoring.md
    for the theoretical framework (Bissiri, Holmes, & Walker 2016) and
    ../../../reading-literature/pymc-tutorial/AR6_FOURIER_POOLED_RESULTS.md
    for the experimental results motivating this choice.

Pooling toggle, selected via `theta_pooling`/`fourier_pooling` (each
"shared" or "none"; `sigma` is always per-location, matching production):
  - "shared": one set of AR(6)/Fourier coefficients for every location.
  - "none": each location gets its own independent coefficients.
"""

from __future__ import annotations

# --- Local macOS dev-machine fix: PyTensor auto-adds a "-ld64" linker flag on
# macOS 15+ to work around an old Xcode-15 linker bug; on macOS 26 with
# current Command Line Tools that flag no longer exists and breaks C
# compilation outright ("ld: library 'd64' not found"). This patch strips it
# and is a harmless no-op on any other platform (including Unity's Linux
# environment, where this code path is never reached). See
# reading-literature/pymc-tutorial/pymc_tutorial.py's first cell for the
# original diagnosis of this issue.
import sys as _sys

if _sys.platform == "darwin":
    from pytensor.link.c.cmodule import GCC_compiler as _GCC_compiler

    _orig_compile_args = _GCC_compiler.compile_args

    def _patched_compile_args(march_flags=True):
        return [f for f in _orig_compile_args(march_flags) if f != "-ld64"]

    _GCC_compiler.compile_args = staticmethod(_patched_compile_args)

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
from scipy import stats

FOURIER_K = 2
AR_ORDER = 6
FLUSIGHT_QUANTILE_LEVELS = np.array(
    [0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
     0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99]
)
_ALPHA_PAIRS = sorted({round(2 * lv, 6) if lv < 0.5 else round(2 * (1 - lv), 6)
                        for lv in FLUSIGHT_QUANTILE_LEVELS if lv != 0.5})
WIS_ALPHA = np.array(_ALPHA_PAIRS)
K_WIS = len(WIS_ALPHA)


# --------------------------------------------------------------------------
# Data transform
# --------------------------------------------------------------------------
#
# WISAR6_fourthroot fits and computes its WIS loss entirely on this
# transformed scale (never on natural counts) -- so the choice of transform
# is also, in effect, a choice of *which* WIS FluSight metric the model is
# aligned with. This version uses `(rate + offset)^0.25`, matching:
#
#   - idmodels.transforms.FourthRootTransform, used by this hub's production
#     AR6_pooled/AR6_fourierP_thetaP models (and this model's own earlier
#     version, before the log1p variant in ../WISAR6/model.py).
#   - the natural-scale FluSight WIS metric as approximated on the
#     fourth-root scale, rather than ../WISAR6/model.py's alignment with the
#     hubverse "log WIS" convention (log(y+1), per Bosse, Abbott, Cori, van
#     Leeuwen, Bracher, & Funk 2023, medRxiv 2023.01.23.23284722).
#
# POWER_TRANSFORM_OFFSET (0.01) matches idmodels.constants.
# POWER_TRANSFORM_OFFSET exactly, so this model's transform is a drop-in
# match for production's, not just the same functional form with a
# different constant.
#
# This is *not* a purely cosmetic choice vs. the log1p sibling: the log and
# fourth-root transforms compress large counts differently (log compresses
# far more aggressively, especially near zero -- see the discussion in
# reading-literature/scoring_estimation/nimble_wis_flusight.md's
# "Natural-Scale WIS vs. Log-Scale WIS" section), so switching transforms
# changes what the fitted AR(6)/Fourier dynamics look like, not just which
# scale the loss is reported on. See
# reading-literature/pymc-tutorial/AR6_ALL_LOCATIONS_RESULTS.md for the
# fourth-root-transform results this version's outputs should be compared
# against, and AR6_ALL_LOCATIONS_LOGWIS_RESULTS.md for the log1p-sibling's.

POWER_TRANSFORM_OFFSET = 0.01


def fourth_root_transform(values: np.ndarray, population: float,
                           offset: float = POWER_TRANSFORM_OFFSET) -> np.ndarray:
    """rate-per-100k, then (rate + offset)^0.25 -- matching
    idmodels.transforms.FourthRootTransform, applied here at fit time rather
    than only at evaluation time."""
    rate = values / population * 100_000
    return (rate + offset) ** 0.25


def fit_scale_params(r_train: np.ndarray) -> tuple[float, float]:
    """Returns (q95, mean_of_scaled), computed from training data only."""
    q95 = np.quantile(r_train, 0.95)
    mean_scaled = np.mean(r_train / q95)
    return q95, mean_scaled


def apply_center_scale(r: np.ndarray, q95: float, mean_scaled: float) -> np.ndarray:
    return r / q95 - mean_scaled


def invert_center_scale(z: np.ndarray, q95: float, mean_scaled: float) -> np.ndarray:
    return (z + mean_scaled) * q95


def invert_fourth_root_transform(r: np.ndarray, population: float,
                                  offset: float = POWER_TRANSFORM_OFFSET) -> np.ndarray:
    rate = np.clip(r, 0, None) ** 4 - offset
    return np.maximum(rate, 0.0) / 100_000 * population


def fourier_design_from_dates(dates) -> np.ndarray:
    """Fourier design matrix on calendar day-of-year, period 365.25 days."""
    doy = pd.DatetimeIndex(dates).dayofyear.values.astype(float)
    cols = []
    for k in range(1, FOURIER_K + 1):
        cols.append(np.sin(2 * np.pi * k * doy / 365.25))
        cols.append(np.cos(2 * np.pi * k * doy / 365.25))
    return np.column_stack(cols)


# --------------------------------------------------------------------------
# WIS loss (used only when objective="wis")
# --------------------------------------------------------------------------

def wis_from_quantiles(y, median_q, lower, upper, alpha, K):
    total = 0.5 * pt.abs(y - median_q)
    denom = 0.5
    for k in range(K):
        w_k = alpha[k] / 2
        pen = upper[:, k] - lower[:, k]
        pen = pen + (2 / alpha[k]) * pt.maximum(lower[:, k] - y, 0)
        pen = pen + (2 / alpha[k]) * pt.maximum(y - upper[:, k], 0)
        total = total + w_k * pen
        denom = denom + w_k
    return total / denom


def gaussian_wis_loss(mu, sigma_vec, y, alpha, K):
    """sigma_vec: per-observation sigma (already indexed by each obs's location)."""
    tau_lower = alpha / 2
    tau_upper = 1 - alpha / 2
    z_lower = pt.sqrt(2) * pt.erfinv(2 * tau_lower - 1)
    z_upper = pt.sqrt(2) * pt.erfinv(2 * tau_upper - 1)
    lower = mu[:, None] + sigma_vec[:, None] * z_lower[None, :]
    upper = mu[:, None] + sigma_vec[:, None] * z_upper[None, :]
    return pt.sum(wis_from_quantiles(y, mu, lower, upper, alpha, K))


def studentt_wis_loss(mu, sigma_vec, nu, y, alpha, K):
    """Student-t analogue of gaussian_wis_loss: quantiles obtained via
    pm.icdf on a StudentT distribution rather than the closed-form Gaussian
    erfinv formula. `nu` (degrees of freedom) MUST be a fixed Python/numpy
    float or a non-sampled pytensor constant here, not a free RV -- PyTensor's
    betaincinv (which pm.icdf's StudentT implementation uses internally) has
    no implemented gradient with respect to its degrees-of-freedom argument
    (confirmed directly: NullTypeGradError on pt.grad w.r.t. nu through
    pm.icdf(StudentT.dist(...))), so NUTS cannot sample nu through this path.
    Gradients w.r.t. mu and sigma_vec *are* implemented and were confirmed
    correct against scipy.stats.t.ppf. See fit_model's docstring for how nu
    is set."""
    tau_lower = alpha / 2
    tau_upper = 1 - alpha / 2
    dist = pm.StudentT.dist(nu=nu, mu=0.0, sigma=1.0)
    t_lower = pm.icdf(dist, tau_lower)
    t_upper = pm.icdf(dist, tau_upper)
    lower = mu[:, None] + sigma_vec[:, None] * t_lower[None, :]
    upper = mu[:, None] + sigma_vec[:, None] * t_upper[None, :]
    return pt.sum(wis_from_quantiles(y, mu, lower, upper, alpha, K))


def wis_reference_from_dict(y, quantiles_dict, levels):
    """Plain-numpy WIS from a {level: value} dict, matching Bracher et al. 2021."""
    m = quantiles_dict[0.5]
    total = 0.5 * abs(y - m)
    denom = 0.5
    alphas = sorted({round(2 * lv, 6) if lv < 0.5 else round(2 * (1 - lv), 6)
                      for lv in levels if lv != 0.5})
    for alpha in alphas:
        lo = quantiles_dict[round(alpha / 2, 6)]
        up = quantiles_dict[round(1 - alpha / 2, 6)]
        w = alpha / 2
        pen = up - lo
        if y < lo:
            pen += (2 / alpha) * (lo - y)
        if y > up:
            pen += (2 / alpha) * (y - up)
        total += w * pen
        denom += w
    return total / denom


def calibrate_loss_scale_pooled(y, loc_idx, F_target, F_lags, z_lags,
                                 target_relative_scale: float = 2.0) -> float:
    """Deterministic OLS-based heuristic for the WIS loss-scale hyperparameter w,
    so exp(-w * WIS) sits on a comparable footing to a Gaussian log-likelihood.
    See generalized_bayesian_scoring.md's "Gibbs Posteriors" section for why
    this scale has no automatically-correct value and remains an open problem;
    this is a rough calibration, not a solved one.

    Always uses gaussian_wis_loss here regardless of this model's own
    innovation_dist setting -- this is only a rough order-of-magnitude
    calibration at a cheap OLS plug-in estimate, not the actual fitted
    model, so the extra complexity of a Student-t version isn't warranted."""
    n_loc = loc_idx.max() + 1
    loc_dummies = np.eye(n_loc)[loc_idx]
    F_lags_flat = F_lags.reshape(F_lags.shape[0], -1)
    X = np.column_stack([loc_dummies, F_target, F_lags_flat, z_lags])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    mu_hat = X @ beta
    resid = y - mu_hat
    sigma_hat_by_loc = np.clip(
        np.array([resid[loc_idx == l].std() for l in range(n_loc)]), 1e-3, None
    )
    sigma_hat_obs = sigma_hat_by_loc[loc_idx]
    loglik = stats.norm.logpdf(y, mu_hat, sigma_hat_obs).sum()
    with pm.Model():
        loss_val = gaussian_wis_loss(
            pt.constant(mu_hat), pt.constant(sigma_hat_obs), pt.constant(y),
            pt.constant(WIS_ALPHA), K_WIS
        ).eval()
    return float(target_relative_scale * abs(loglik) / loss_val)


# --------------------------------------------------------------------------
# Design matrix construction
# --------------------------------------------------------------------------

def build_full_design(z_by_loc: dict, dates_by_loc: dict, locs: list, p: int = AR_ORDER):
    """Stack all locations' training series into flat arrays for a pooled/unpooled fit.
    Returns (y, loc_idx, F_target, F_lags, z_lags)."""
    y_all, loc_idx_all = [], []
    F_target_all, F_lags_all, z_lags_all = [], [], []
    for li, loc in enumerate(locs):
        z = z_by_loc[loc]
        F = fourier_design_from_dates(dates_by_loc[loc])
        T = len(z)
        n = T - p
        y_all.append(z[p:])
        loc_idx_all.append(np.full(n, li))
        F_target_all.append(F[p:])
        F_lags = np.stack([F[p - j: T - j] for j in range(1, p + 1)], axis=1)
        z_lags = np.stack([z[p - j: T - j] for j in range(1, p + 1)], axis=1)
        F_lags_all.append(F_lags)
        z_lags_all.append(z_lags)
    return (np.concatenate(y_all), np.concatenate(loc_idx_all),
            np.concatenate(F_target_all, axis=0), np.concatenate(F_lags_all, axis=0),
            np.concatenate(z_lags_all, axis=0))


def build_pooled_mean(phi, fourier_coefs, mu0, theta_pooling, fourier_pooling,
                       loc_idx, F_target, F_lags, z_lags):
    mu0_obs = mu0[loc_idx]
    if fourier_pooling == "shared":
        S_target = pt.dot(F_target, fourier_coefs)
        S_lags = pt.dot(
            F_lags.reshape((-1, F_lags.shape[-1])), fourier_coefs
        ).reshape(F_lags.shape[:-1])
    else:
        fc_obs = fourier_coefs[loc_idx]
        S_target = pt.sum(F_target * fc_obs, axis=1)
        S_lags = pt.sum(F_lags * fc_obs[:, None, :], axis=2)

    dev_lags = z_lags - mu0_obs[:, None] - S_lags

    if theta_pooling == "shared":
        ar_term = pt.dot(dev_lags, phi)
    else:
        phi_obs = phi[loc_idx]
        ar_term = pt.sum(dev_lags * phi_obs, axis=1)

    return mu0_obs + S_target + ar_term


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------

def fit_model(
    y, loc_idx, F_target, F_lags, z_lags, n_loc: int,
    objective: str = "wis",
    theta_pooling: str = "shared",
    fourier_pooling: str = "shared",
    loss_scale: float = 1.0,
    draws: int = 2000,
    tune: int = 2000,
    chains: int = 2,
    cores: int = 1,
    seed: int = 0,
    target_accept: float = 0.99,
    sigma_prior_scale: float = 1.0,
    theta_sd_prior_scale: float = 1.0,
    fourier_beta_sd_prior_scale: float = 1.0,
    innovation_dist: str = "studentt",
    t_dof: float = 5.0,
    callback=None,
):
    """Fit the pooled AR(6)+Fourier model.

    NOTE on `cores=1`: PyMC's default multiprocessing (one OS subprocess per
    chain) was found to crash with EOFError at 53-location scale during
    development (confirmed via isolated reproduction: cores=1 succeeds
    immediately with the identical model; the crash is a multiprocessing
    quirk on the development machine, not a memory or model-correctness
    issue). `cores=1` runs chains sequentially in a single process --
    slower than parallel chains, but avoids the crash. This default was
    chosen on macOS, where this specific EOFError is a known class of issue
    with Python multiprocessing + compiled C extensions; it has not been
    observed to occur on Linux, which uses a different (fork-based) process
    start method and doesn't share this failure mode. On Unity (Linux),
    main.py sets MCMC_CORES=2 via the WISAR6_MCMC_CORES env var, which
    should roughly halve sampling wall-clock time (2 chains run in parallel
    instead of sequentially) -- but is not yet validated at full production
    scale on Unity itself, so the first Unity run should be watched for the
    same crash before trusting it unattended across a full SLURM array.

    Priors: mu0 ~ Normal(0, 0.1^2), unchanged from the WISAR6/WISAR6_fourthroot
    siblings. phi, sigma, and the Fourier coefficients each instead now use a
    learned (data-adaptive) prior scale, mirroring production SARIX's own
    prior structure exactly (confirmed against the installed `sarix` package
    source, `sarix.py`'s `model()` method):

        sigma         ~ HalfCauchy(sigma_prior_scale)   -- directly the
                        innovation scale, same as production (not a
                        hyperprior on a further scale parameter)
        theta_sd      ~ HalfCauchy(theta_sd_prior_scale);  phi ~ Normal(0, theta_sd)
        fourier_beta_sd ~ HalfCauchy(fourier_beta_sd_prior_scale);
                        fourier_coefs ~ Normal(0, fourier_beta_sd)

    This replaces the WISAR6/WISAR6_fourthroot siblings' fixed
    phi ~ Normal(0, 0.15^2), sigma ~ HalfNormal(0, 0.5^2), and
    fourier_coefs ~ Normal(0, 1^2) -- each of which was a fixed constant with
    no hyperprior. All three *_prior_scale hyperparameters default to 1.0,
    matching SARIX's own defaults exactly (sarix.SARIX.__init__'s
    sigma_prior_scale, theta_sd_prior_scale, fourier_beta_sd_prior_scale).

    This change was motivated by a WIS evaluation of WISAR6_fourthroot
    against AR6_pooled/AR6_fourierP_thetaP that found elevated
    underprediction both near a season's peak/shoulders (consistent with the
    Fourier term being too rigid to track the real curve's acceleration) and,
    unexpectedly, during calm off-season weeks too (not explained by
    seasonal-curve rigidity alone, and more consistent with sigma/phi's fixed
    priors being systematically too tight even when nothing seasonal is
    happening) -- see
    ../../model-output/UMass-WISAR6_fourthroot/PERFORMANCE_NOTES.md for the
    full diagnosis and ../../peak_proximity_analysis.R for the supporting
    analysis.

    The *fixed* phi prior in the WISAR6/WISAR6_fourthroot siblings was
    originally tightened specifically to resolve a near-unit-root divergence
    issue in the US series (its AR(6) dynamics land at a characteristic root
    magnitude of ~0.94, right at the edge of stationarity, once the seasonal
    component absorbs most of the systematic variation -- see
    ar6_fourier_pooled_flusight.py and AR6_FOURIER_POOLED_RESULTS.md for the
    full diagnosis). Loosening phi's prior back to a learned scale risks
    reintroducing that instability; `target_accept=0.99` (unchanged from the
    siblings) is the first line of defense, and if divergences reappear in
    practice the next lever is tightening `theta_sd_prior_scale` below its
    1.0 default (a weakly-informative upper bound on the hyperprior) rather
    than reintroducing a fixed, hand-picked phi scale.

    `innovation_dist` selects the innovation distribution eps_{ell,t}:
      - "gaussian": eps ~ Normal(0, sigma_ell^2), matching every prior
        WISAR6 sibling and production SARIX.
      - "studentt" (default here): eps ~ sigma_ell * StudentT(nu=t_dof),
        a heavier-tailed alternative intended to better accommodate the
        sharp, occasionally extreme week-to-week swings a linear Gaussian
        AR(6) cannot capture -- the same robustness motivation originally
        cited for the WIS-loss objective itself (see
        AR6_ALL_LOCATIONS_RESULTS.md section 1.1), applied here to the
        *distributional* assumption instead. `t_dof` MUST be a fixed
        constant, not sampled, when objective="wis": PyTensor's
        betaincinv (which the Student-t WIS-quantile construction in
        studentt_wis_loss uses via pm.icdf) has no implemented gradient
        with respect to the degrees-of-freedom argument (confirmed
        directly -- see studentt_wis_loss's docstring), so NUTS cannot
        sample nu through that path; t_dof=5 is a commonly-used default
        for moderately heavy tails.
    """
    loc_idx_pt = pt.constant(loc_idx)
    F_target_pt = pt.constant(F_target)
    F_lags_pt = pt.constant(F_lags)
    z_lags_pt = pt.constant(z_lags)
    y_pt = pt.constant(y)

    with pm.Model() as model:
        mu0 = pm.Normal("mu0", 0, 0.1, shape=n_loc)
        sigma = pm.HalfCauchy("sigma", sigma_prior_scale, shape=n_loc)

        theta_sd = pm.HalfCauchy("theta_sd", theta_sd_prior_scale)
        if theta_pooling == "shared":
            phi = pm.Normal("phi", 0, theta_sd, shape=AR_ORDER)
        else:
            phi = pm.Normal("phi", 0, theta_sd, shape=(n_loc, AR_ORDER))

        fourier_beta_sd = pm.HalfCauchy("fourier_beta_sd", fourier_beta_sd_prior_scale)
        if fourier_pooling == "shared":
            fourier_coefs = pm.Normal("fourier_coefs", 0, fourier_beta_sd, shape=2 * FOURIER_K)
        else:
            fourier_coefs = pm.Normal("fourier_coefs", 0, fourier_beta_sd, shape=(n_loc, 2 * FOURIER_K))

        mu = build_pooled_mean(phi, fourier_coefs, mu0, theta_pooling, fourier_pooling,
                                loc_idx_pt, F_target_pt, F_lags_pt, z_lags_pt)
        sigma_obs = sigma[loc_idx_pt]

        if innovation_dist not in ("gaussian", "studentt"):
            raise ValueError(f"Unknown innovation_dist: {innovation_dist!r}")

        if objective == "likelihood":
            if innovation_dist == "gaussian":
                pm.Normal("y_obs", mu, sigma_obs, observed=y_pt)
            else:
                pm.StudentT("y_obs", nu=t_dof, mu=mu, sigma=sigma_obs, observed=y_pt)
        elif objective == "wis":
            if innovation_dist == "gaussian":
                loss = gaussian_wis_loss(mu, sigma_obs, y_pt, pt.constant(WIS_ALPHA), K_WIS)
            else:
                loss = studentt_wis_loss(mu, sigma_obs, t_dof, y_pt, pt.constant(WIS_ALPHA), K_WIS)
            pm.Potential("wis_potential", -loss_scale * loss)
        else:
            raise ValueError(f"Unknown objective: {objective!r}")

        idata = pm.sample(draws, tune=tune, chains=chains, cores=cores,
                           target_accept=target_accept, random_seed=seed,
                           progressbar=False, callback=callback)
    return model, idata


# --------------------------------------------------------------------------
# Forecasting (Monte Carlo simulation, matching how a pooled model with
# shared parameters must generate forecasts -- no single closed form once
# phi is shared across locations with different mu0/sigma)
# --------------------------------------------------------------------------

def simulate_forecast_draws(idata, pooling: str, loc_index: int, z_history: np.ndarray,
                             dates_history, h_max: int, n_param_draws: int = 200,
                             n_innov_per_param: int = 5, seed: int = 0,
                             innovation_dist: str = "studentt", t_dof: float = 5.0) -> np.ndarray:
    """Returns array of shape (n_param_draws * n_innov_per_param, h_max).

    `innovation_dist`/`t_dof` must match the values used in the `fit_model`
    call that produced `idata`, so the simulated forecast distribution
    matches the fitted one."""
    if innovation_dist not in ("gaussian", "studentt"):
        raise ValueError(f"Unknown innovation_dist: {innovation_dist!r}")
    rng = np.random.default_rng(seed)
    post = idata.posterior
    n_chains, n_draws = post.sizes["chain"], post.sizes["draw"]
    flat_idx = rng.choice(n_chains * n_draws, size=n_param_draws, replace=False)
    chain_idx, draw_idx = np.unravel_index(flat_idx, (n_chains, n_draws))

    mu0_all = post["mu0"].values[chain_idx, draw_idx, loc_index]
    sigma_all = post["sigma"].values[chain_idx, draw_idx, loc_index]
    if pooling == "shared":
        phi_all = post["phi"].values[chain_idx, draw_idx, :]
        fc_all = post["fourier_coefs"].values[chain_idx, draw_idx, :]
    else:
        phi_all = post["phi"].values[chain_idx, draw_idx, loc_index, :]
        fc_all = post["fourier_coefs"].values[chain_idx, draw_idx, loc_index, :]

    last_dates = pd.DatetimeIndex(dates_history[-AR_ORDER:])
    future_dates = pd.date_range(dates_history[-1], periods=h_max + 1, freq="W-SAT")[1:]
    F_future = fourier_design_from_dates(future_dates)
    F_recent = fourier_design_from_dates(last_dates)

    all_sims = np.zeros((n_param_draws, n_innov_per_param, h_max))
    for i in range(n_param_draws):
        mu0, sigma, phi, fc = mu0_all[i], sigma_all[i], phi_all[i], fc_all[i]
        S_recent = F_recent @ fc
        dev_recent = z_history[-AR_ORDER:] - mu0 - S_recent
        for j in range(n_innov_per_param):
            dev = list(dev_recent)
            for h in range(h_max):
                ar_mean = np.dot(phi, np.array(dev[-AR_ORDER:][::-1]))
                S_h = F_future[h] @ fc
                if innovation_dist == "gaussian":
                    innovation = rng.normal(0, sigma)
                else:
                    innovation = rng.standard_t(t_dof) * sigma
                z_h = mu0 + S_h + ar_mean + innovation
                all_sims[i, j, h] = z_h
                dev.append(z_h - mu0 - S_h)
    return all_sims.reshape(-1, h_max)
