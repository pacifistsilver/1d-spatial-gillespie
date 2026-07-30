<<<<<<< HEAD
<<<<<<< HEAD
"""Bayesian models for the monomer and heterodimer promoters.

These were ABC-SMC models: a Gillespie simulator compared to the data through
``summary_stat``, which collapsed both to ``[mean, fano]`` and accepted anything
within an epsilon ball. They now use the exact stationary likelihood from
:mod:`stochtf.inference.likelihood`, which scores the whole distribution
deterministically. The sampler is still SMC, but it is no longer *approximate*
Bayesian computation -- it targets the true posterior.

Two modelling points that the ABC version got wrong:

Identifiability
    Stationary counts determine only the rates *relative to* the mRNA
    degradation rate: the distribution depends on
    (alpha_s/gamma, beta_s/gamma, alpha_n/gamma, beta_n/gamma, k_y/gamma) and
    not on gamma separately. gamma is therefore fixed at 1, which sets the time
    unit to one mRNA lifetime, and every other rate is inferred in those units.

Transcription rate
    The old simulators took ``k_y=0.01, gamma_y=0.005`` as *default arguments*
    that ``pm.Simulator`` never overrode, so k_y/gamma was pinned at 2 and the
    model mean could not exceed 2 (OR) or 4 (ADD). The observed means run from
    29.7 (sox2) to 251.5 (esrrb), so no parameter setting could fit the data.
    k_y is now inferred.
=======
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
"""ABC-SMC model definitions for the monomer and heterodimer promoters.

Both models previously duplicated ``_data_setter``,
``_generate_and_preprocess_model_data``, the two config staticmethods, ``fit``
and a 37-line ``save`` verbatim; the shared behaviour now lives on
:class:`_SMCModel` and each subclass supplies only its priors and simulator.
<<<<<<< HEAD
>>>>>>> 96a2b5c (refactor: restructure into an installable package for publication)
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
"""

import json
from typing import Dict

<<<<<<< HEAD
<<<<<<< HEAD
import numpy as np
import pymc as pm
import pytensor.tensor as pt
from pymc_extras.model_builder import ModelBuilder
from pytensor.graph.basic import Apply
from pytensor.graph.op import Op

from stochtf.inference.likelihood import log_likelihood, prepare_counts

#: Time unit: rates are expressed relative to the mRNA degradation rate, which
#: stationary count data cannot identify separately.
GAMMA = 1.0

#: Log-normal prior scales, in units of gamma. Broad enough to cover the
#: observed dynamic range (means from ~30 to ~250 molecules).
DEFAULT_MODEL_CONFIG: Dict = {
    "alpha_s_mu": 4, "alpha_s_sigma": 1,
    "alpha_n_mu": 4, "alpha_n_sigma": 1,
    "beta_s_mu": 2, "beta_s_sigma": 1,
    "beta_n_mu": 3, "beta_n_sigma": 1,
    "k_y_mu": 3, "k_y_sigma": 1.0,
}


DEFAULT_SAMPLER_CONFIG: Dict = {
    "draws": 1000,
    "chains": 4,
}


class StationaryLogLike(Op):
    """Gradient-free PyTensor Op wrapping the exact stationary log-likelihood.

    The likelihood needs a linear solve per parameter proposal, so there is no
    tractable gradient. That is fine: ``pm.sample_smc`` uses gradient-free
    Metropolis kernels within each stage.

    This is a hand-written Op rather than ``pytensor.wrap_py`` (formerly
    ``as_op``) because the Op has to survive pickling. ``wrap_py`` reduces by
    *name lookup* -- it pickles as ``getattr(module, func.__name__)`` -- so an
    Op built inside a factory, closing over the observed counts, cannot be
    unpickled: the closure has no module-level name. ``pm.sample_smc`` runs
    chains in separate processes and cloudpickles the kernel, so that failed
    with "Can't pickle wrap_py(), not found as stochtf.inference.models._logp".

    A normal Op subclass pickles by state instead: the class resolves by
    reference and ``counts``/``model_name`` travel by value.

    No ``__props__`` is declared, so equality falls back to object identity.
    That is deliberate -- two instances holding different datasets must never
    compare equal, or PyTensor's graph merging could substitute one for the
    other.
    """

    def __init__(self, counts, model_name):
        self.counts = prepare_counts(counts)
        self.model_name = model_name

    def make_node(self, alpha_s, alpha_n, beta_s, beta_n, k_y):
        inputs = [pt.as_tensor_variable(v)
                  for v in (alpha_s, alpha_n, beta_s, beta_n, k_y)]
        return Apply(self, inputs, [pt.scalar(dtype="float64")])

    def perform(self, node, inputs, output_storage):
        alpha_s, alpha_n, beta_s, beta_n, k_y = (float(v) for v in inputs)
        value = log_likelihood(
            self.counts,
            alpha_s, beta_s, alpha_n, beta_n, k_y, GAMMA,
            model=self.model_name,
        )
        output_storage[0][0] = np.asarray(value, dtype="float64")


class _ExactModel(ModelBuilder):
    """Shared plumbing. Subclasses set :attr:`gate_model` only."""

    #: Key into stochtf.analytical.pgf.MODEL_GATE.
    gate_model = None
    version = "1.0"

    def _data_setter(self, X, y=None):
        pass  # ModelBuilder requires this method; X/y are not standard inputs
=======
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
import polars as pl
import pymc as pm
from pymc_extras.model_builder import ModelBuilder

from stochtf.inference.abc_smc import (
    pymc_fast_simulator_het,
    pymc_fast_simulator_monomer,
    summary_stat,
)

#: Prior scales shared by both models.
DEFAULT_MODEL_CONFIG: Dict = {
    "alpha_s_sigma": 1.0,
    "alpha_n_sigma": 1.0,
    "beta_s_sigma": 0.1,
    "beta_n_sigma": 0.3,
    "gamma_y_sigma": 0.01,
    "k_y_sigma": 0.5,
}

DEFAULT_SAMPLER_CONFIG: Dict = {
    "draws": 500,
    "tune": 1000,
    "chains": 5,
    "target_accept": 0.95,
}

#: Tolerance for the ABC acceptance kernel.
EPSILON = 2.0


class _SMCModel(ModelBuilder):
    """Shared ABC-SMC plumbing. Subclasses define ``_prior`` and ``_simulator``."""

    #: pm distribution used for all four rate priors.
    _prior = None
    #: numba-compiled simulator passed to pm.Simulator.
    _simulator = None

    def _data_setter(self, X, y=None):
        pass  # ModelBuilder requires this method; X/y are not used as standard inputs
<<<<<<< HEAD
>>>>>>> 96a2b5c (refactor: restructure into an installable package for publication)
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)

    def _generate_and_preprocess_model_data(self, X, y=None):
        pass  # ModelBuilder requires this method

    @staticmethod
    def get_default_model_config() -> Dict:
<<<<<<< HEAD
<<<<<<< HEAD
=======
        """Prior scales used when no model_config is supplied at construction."""
>>>>>>> 96a2b5c (refactor: restructure into an installable package for publication)
=======
        """Prior scales used when no model_config is supplied at construction."""
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
        return dict(DEFAULT_MODEL_CONFIG)

    @staticmethod
    def get_default_sampler_config() -> Dict:
<<<<<<< HEAD
<<<<<<< HEAD
        return dict(DEFAULT_SAMPLER_CONFIG)

    def build_model(self, X=None, y=None, **kwargs):
        if y is None:
            raise ValueError("observed counts are required to build the model")
        cfg = self.model_config
        logp_op = StationaryLogLike(y, self.gate_model)

        with pm.Model() as self.model:
            alpha_s = pm.LogNormal("alpha_s", mu=cfg["alpha_s_mu"],
                                   sigma=cfg["alpha_s_sigma"])
            alpha_n = pm.LogNormal("alpha_n", mu=cfg["alpha_n_mu"],
                                   sigma=cfg["alpha_n_sigma"])
            beta_s = pm.LogNormal("beta_s", mu=cfg["beta_s_mu"],
                                  sigma=cfg["beta_s_sigma"])
            beta_n = pm.LogNormal("beta_n", mu=cfg["beta_n_mu"],
                                  sigma=cfg["beta_n_sigma"])
            k_y = pm.LogNormal("k_y", mu=cfg["k_y_mu"], sigma=cfg["k_y_sigma"])

            # Grouped with the data log-probability, so SMC tempers it correctly.
            pm.Potential("likelihood",
                         logp_op(alpha_s, alpha_n, beta_s, beta_n, k_y))

        return self.model

    def fit(self, data, sampler_config: dict = None, **kwargs):
        """Fit by Sequential Monte Carlo over the exact posterior."""
        if sampler_config is None:
            sampler_config = self.sampler_config

        self.build_model(y=data)

        with self.model:
            self.idata = pm.sample_smc(
                draws=sampler_config["draws"],
                chains=sampler_config.get("chains", 4),
                **kwargs,
            )
            prior = pm.sample_prior_predictive(draws=500, random_seed=500)
            if "prior" in prior:
                self.idata["prior"] = prior["prior"]

=======
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
        """Sampler settings used when no sampler_config is supplied."""
        return dict(DEFAULT_SAMPLER_CONFIG)

    def build_model(self, X=None, y=None, **kwargs):
        with pm.Model() as self.model:
            cfg = self.model_config
            alpha_s = self._prior("alpha_s", sigma=cfg["alpha_s_sigma"])
            alpha_n = self._prior("alpha_n", sigma=cfg["alpha_n_sigma"])
            beta_s = self._prior("beta_s", sigma=cfg["beta_s_sigma"])
            beta_n = self._prior("beta_n", sigma=cfg["beta_n_sigma"])

            pm.Simulator(
                "sim",
                self._simulator,
                params=(alpha_s, alpha_n, beta_s, beta_n),
                sum_stat=summary_stat,
                epsilon=EPSILON,
                observed=y,
            )

    def fit(self, data: pl.DataFrame, sampler_config: dict = None, **kwargs):
        """Fit by Sequential Monte Carlo instead of ModelBuilder's default NUTS."""
        if sampler_config is None:
            sampler_config = self.sampler_config

        self.build_model(self.model_config, y=data)

        with self.model:
            self.idata = pm.sample_smc(draws=sampler_config["draws"], **kwargs)

            prior_distribution = pm.sample_prior_predictive(draws=500, random_seed=500)
            log_likelihood = pm.compute_log_likelihood(self.idata)
            log_prior = pm.stats.compute_log_prior(self.idata)

            if "prior" in prior_distribution:
                self.idata["prior"] = prior_distribution["prior"]
            if "prior_predictive" in prior_distribution:
                self.idata["prior_predictive"] = prior_distribution["prior_predictive"]
            self.idata["log_likelihood"] = log_likelihood["log_likelihood"]
            self.idata["log_prior"] = log_prior["log_prior"]
<<<<<<< HEAD
>>>>>>> 96a2b5c (refactor: restructure into an installable package for publication)
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
        return self.idata

    def save(self, fname: str):
        """Sanitise the DataTree before delegating to ModelBuilder's saver.

        NetCDF cannot represent object-dtype arrays, and ModelBuilder.load()
        raises KeyError if its metadata attrs are absent, so both are fixed up
        here first.
        """
        if getattr(self, "idata", None) is not None:
            attrs = self.idata.attrs
            if "model_config" not in attrs:
                attrs["model_config"] = json.dumps(self.model_config)
            if "sampler_config" not in attrs:
                attrs["sampler_config"] = json.dumps(
                    self.sampler_config
                    if hasattr(self, "sampler_config")
                    else self.get_default_sampler_config()
                )
            if "version" not in attrs:
                attrs["version"] = self.version
            if "model_type" not in attrs:
                attrs["model_type"] = self.model_type

            for group_name in self.idata.groups:
                group = self.idata[group_name]
                for var_name in list(group.data_vars.keys()):
                    var = group[var_name]
                    if var.dtype == object:
                        try:
                            # Force integer/float mixes into pure floats
                            group[var_name] = var.astype(float)
                        except Exception:
                            # Sequence arrays cannot be coerced; drop them
                            del group[var_name]

        super().save(fname)


<<<<<<< HEAD
<<<<<<< HEAD
class MonomerModel(_ExactModel):
    """Independent SOX2 and NANOG sites, additive transcription rate.

    Only k_y and one occupancy-like combination are identifiable from
    stationary counts; see :class:`TelegraphModel` and
    :mod:`stochtf.inference.identifiability`.
    """

    model_type = "MonomerModel"
    gate_model = "monomer"


class DimerModel(_ExactModel):
    """Heterodimer occupancy: any bound site gives the full transcription rate.

    Same caveat as :class:`MonomerModel`: the four switching rates are not
    determined by stationary counts.
    """

    model_type = "DimerModel"
    gate_model = "dimer"


#: Site n is switched off, reducing the promoter to a single site.
SITE_OFF = 1e-9

#: Priors for the single-site model, in units of gamma.
TELEGRAPH_MODEL_CONFIG: Dict = {
    "alpha_s_mu": 0.0, "alpha_s_sigma": 2.0,
    "beta_s_mu": 0.0, "beta_s_sigma": 2.0,
    "k_y_mu": 3.0, "k_y_sigma": 2.0,
}


class TelegraphModel(_ExactModel):
    """Single-site (telegraph) promoter: alpha, beta and k_y only.

    This is the one variant whose parameters stationary counts can actually
    determine -- around 7%, 10% and 1% at n = 800, matching published
    Poisson-Beta fits. It is therefore the model to use for an end-to-end
    recovery test of the likelihood and sampler.

    Site n is pinned off rather than given its own free rates, which both
    reduces the promoter to two states and removes the site-exchange symmetry.
    """

    model_type = "TelegraphModel"
    gate_model = "dimer"

    @staticmethod
    def get_default_model_config() -> Dict:
        return dict(TELEGRAPH_MODEL_CONFIG)

    def build_model(self, X=None, y=None, **kwargs):
        if y is None:
            raise ValueError("observed counts are required to build the model")
        cfg = self.model_config
        logp_op = StationaryLogLike(y, self.gate_model)

        with pm.Model() as self.model:
            alpha_s = pm.LogNormal("alpha_s", mu=cfg["alpha_s_mu"],
                                   sigma=cfg["alpha_s_sigma"])
            beta_s = pm.LogNormal("beta_s", mu=cfg["beta_s_mu"],
                                  sigma=cfg["beta_s_sigma"])
            k_y = pm.LogNormal("k_y", mu=cfg["k_y_mu"], sigma=cfg["k_y_sigma"])

            off = pt.constant(SITE_OFF)
            one = pt.constant(1.0)
            pm.Potential("likelihood", logp_op(alpha_s, off, beta_s, one, k_y))

        return self.model


MODELS = {"monomer": MonomerModel, "dimer": DimerModel,
          "telegraph": TelegraphModel}
=======
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
class MonomerModel(_SMCModel):
    """Independent SOX2 and NANOG binding, half-normal priors."""

    model_type = "MonomerModel"
    version = "0.1"
    _prior = staticmethod(pm.HalfNormal)
    _simulator = staticmethod(pymc_fast_simulator_monomer)


class DimerModel(_SMCModel):
    """Heterodimer binding, log-normal priors."""

    model_type = "DimerModel"
    version = "0.1"
    _prior = staticmethod(pm.LogNormal)
    _simulator = staticmethod(pymc_fast_simulator_het)


MODELS = {"monomer": MonomerModel, "dimer": DimerModel}
<<<<<<< HEAD
>>>>>>> 96a2b5c (refactor: restructure into an installable package for publication)
=======
>>>>>>> 3915fd7 (refactor: restructure into an installable package for publication)
