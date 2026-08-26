"""Which parameters stationary counts can actually determine.

Fitting a model that the data cannot constrain produces posteriors that look
like answers but are really priors. This module measures that directly, before
any sampling, from the Fisher information of the exact stationary distribution.

For iid counts the information about theta is

    I_jk = n * sum_y  (dP(y)/dtheta_j)(dP(y)/dtheta_k) / P(y),

computed here in *log* parameter space, so the resulting standard errors read as
fractional errors: 0.1 means "determined to about 10%".

Result for the two-site promoter:
  A single-site telegraph promoter (alpha, beta, k_y) is well identified --
  around 7%, 10% and 1% at n = 800, matching published Poisson-Beta fits to
  single-cell data. Adding the second site in the same bursty regime sends the
  four switching rates to standard errors of 600-5000%: only k_y and one
  occupancy-like combination survive. Of five log-parameter directions, two are
  determined, one is weak and two are flat. This is a property of the model and
  the data type, not of the sampler or the priors -- stationary counts are a
  one-dimensional summary of a four-state promoter, and no parameter regime
  repairs it.

There is also an exact discrete degeneracy: the likelihood is invariant under
exchanging the two sites, (alpha_s, beta_s) <-> (alpha_n, beta_n), for every
gate. The sites are therefore identifiable at best as an unordered pair, and
marginal posterior means average over both labellings unless the sampler is
given an ordering constraint. See :func:`exchange_symmetry_residual`.
"""

import numpy as np

from stochtf.analytical import pgf

PARAM_NAMES = ("alpha_s", "beta_s", "alpha_n", "beta_n", "k_y")

#: Fractional standard error below which a parameter is called determined.
DETERMINED_SE = 0.5


def _pmf_grid(theta, gate, y_max):
    a_s, b_s, a_n, b_n, k_y = theta
    return pgf.stationary_pmf(a_s, b_s, a_n, b_n, k_y, 1.0, gate, y_max=y_max)


def _default_y_max(theta, gate):
    """Returns a count grid bound covering the bulk of the distribution."""
    a_s, b_s, a_n, b_n, k_y = theta
    mean, var, _ = pgf.moments(a_s, b_s, a_n, b_n, k_y, 1.0, gate)
    return int(mean + 14.0 * np.sqrt(max(var, 1.0)) + 30)


def fisher_information(theta, n_obs, gate="OR", free=None, h=1e-5, y_max=None):
    """Computes the Fisher information in log-parameter space.

    Uses first derivatives only. A finite-difference Hessian of the expected
    log-likelihood is the textbook route but loses too much precision here --
    it returns negative eigenvalues on distributions this flat.

    Args:
        theta: (alpha_s, beta_s, alpha_n, beta_n, k_y), with gamma = 1.
        n_obs: Number of observed cells.
        gate: Promoter logic, one of "OR", "AND", "ADD".
        free: Indices of theta to treat as unknown. The rest are held fixed,
          which is how the single-site case is obtained (freeze site n).
        h: Relative step for the finite-difference derivatives.
        y_max: Count grid bound. Chosen from the moments if None.

    Returns:
        The information matrix, square in the number of free parameters.
    """
    theta = np.asarray(theta, dtype=float)
    if free is None:
        free = range(len(theta))
    free = list(free)
    if y_max is None:
        y_max = _default_y_max(theta, gate)

    log_theta = np.log(theta)
    p0 = np.clip(_pmf_grid(theta, gate, y_max), 1e-300, None)

    derivs = np.empty((len(free), p0.size))
    for row, i in enumerate(free):
        step = np.zeros_like(log_theta)
        step[i] = h
        plus = _pmf_grid(np.exp(log_theta + step), gate, y_max)
        minus = _pmf_grid(np.exp(log_theta - step), gate, y_max)
        derivs[row] = (plus - minus) / (2 * h)

    return n_obs * (derivs / p0) @ derivs.T


def standard_errors(theta, n_obs, gate="OR", free=None, **kwargs):
    """Computes fractional standard errors from the Cramer-Rao bound.

    A pseudo-inverse is used because the information matrix is genuinely
    singular in the flat directions. Values there are lower bounds on an
    already enormous error, so treat anything above ~1 as "not determined"
    rather than as a number.

    Args:
        theta: (alpha_s, beta_s, alpha_n, beta_n, k_y), with gamma = 1.
        n_obs: Number of observed cells.
        gate: Promoter logic, one of "OR", "AND", "ADD".
        free: Indices of theta to treat as unknown.
        **kwargs: Forwarded to :func:`fisher_information`.

    Returns:
        One fractional standard error per free parameter.
    """
    info = fisher_information(theta, n_obs, gate, free, **kwargs)
    return np.sqrt(np.diag(np.linalg.pinv(info)))


def identifiable_directions(theta, n_obs, gate="OR", **kwargs):
    """Decomposes the information matrix into determined directions.

    Args:
        theta: (alpha_s, beta_s, alpha_n, beta_n, k_y), with gamma = 1.
        n_obs: Number of observed cells.
        gate: Promoter logic, one of "OR", "AND", "ADD".
        **kwargs: Forwarded to :func:`fisher_information`.

    Returns:
        A tuple (eigenvalues, eigenvectors) sorted from best- to
        worst-determined. Each eigenvector is a combination of log-parameters:
        large eigenvalues are what the data pins down, near-zero ones are the
        flat directions along which the likelihood is indifferent.
    """
    info = fisher_information(theta, n_obs, gate, **kwargs)
    vals, vecs = np.linalg.eigh(info)
    order = np.argsort(vals)[::-1]
    return vals[order], vecs[:, order]


def describe(theta, n_obs, gate="OR", free=None, names=PARAM_NAMES):
    """Formats a human-readable identifiability report for one theta.

    Args:
        theta: (alpha_s, beta_s, alpha_n, beta_n, k_y), with gamma = 1.
        n_obs: Number of observed cells.
        gate: Promoter logic, one of "OR", "AND", "ADD".
        free: Indices of theta to treat as unknown.
        names: Parameter names used in the report.

    Returns:
        The report as a string.
    """
    free = list(range(len(theta))) if free is None else list(free)
    se = standard_errors(theta, n_obs, gate, free)
    mean, var, fano = pgf.moments(*theta[:4], theta[4], 1.0, gate)

    lines = [f"gate={gate}  n={n_obs}  mean={mean:.2f}  Fano={fano:.2f}",
             "fractional standard errors:"]
    for idx, value in zip(free, se):
        verdict = "determined" if value < DETERMINED_SE else "NOT identified"
        lines.append(f"    {names[idx]:8s} {value * 100:10.1f}%   {verdict}")
    return "\n".join(lines)


def exchange_symmetry_residual(theta, gate="OR", y_max=400):
    """Measures how far the likelihood moves when the two sites are swapped.

    Zero to rounding for every gate: the sites are exchangeable, so the
    posterior has two exactly equivalent modes. Any inference reporting alpha_s
    and alpha_n separately needs an ordering constraint to break it.

    Args:
        theta: (alpha_s, beta_s, alpha_n, beta_n, k_y), with gamma = 1.
        gate: Promoter logic, one of "OR", "AND", "ADD".
        y_max: Count grid bound.

    Returns:
        max |P(theta) - P(theta with the two sites swapped)|.
    """
    a_s, b_s, a_n, b_n, k_y = theta
    p1 = pgf.stationary_pmf(a_s, b_s, a_n, b_n, k_y, 1.0, gate, y_max=y_max)
    p2 = pgf.stationary_pmf(a_n, b_n, a_s, b_s, k_y, 1.0, gate, y_max=y_max)
    return float(np.abs(p1 - p2).max())
