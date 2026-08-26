# stochtf — stochastic modelling of NANOG/SOX2 dynamics in gene expression noise

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![DOI](https://img.shields.io/badge/DOI-10.64898%2F2025.12.03.691924-orange.svg)

Code accompanying:

> G. G. Agsu et al., "Protein-protein interactions drive differences in the
> spatiotemporal dynamics of transcription factors NANOG and SOX2 in naïve
> pluripotent cells," Dec. 2025.
> doi: [10.64898/2025.12.03.691924](https://doi.org/10.64898/2025.12.03.691924)

## What this is

In mouse embryonic stem cells, SOX2 and NANOG regulate the same pluripotency
targets. *How* they share a promoter is not settled, and two architectures are
on the table:

- **HD — two independent sites.** Each factor binds its own site, and the gene
  transcribes whenever at least one site is occupied.
- **M — one contested site.** The factors compete for a single site and are
  never co-bound.

These are different Markov chains, and they leave different fingerprints in the
distribution of mRNA counts across single cells. `stochtf` works out those
fingerprints — in closed form where one exists, by exact solution of the
chemical master equation where it does not — and fits both architectures to
allele-resolved single-cell counts to see which the data support.

```mermaid
graph LR
  subgraph HD["HD — two independent sites"]
    direction LR
    A["00 · silent"] -->|a_s| B["10"]
    A -->|a_n| C["01"]
    B -->|a_n| D["11"]
    C -->|a_s| D
    B -->|b_s| A
    C -->|b_n| A
    D -->|b_n| B
    D -->|b_s| C
  end
  subgraph M["M — one contested site"]
    direction LR
    S["S · SOX2 bound"] -->|b_s| Z["0 · silent"]
    Z -->|a_s| S
    Z -->|a_n| N["N · NANOG bound"]
    N -->|b_n| Z
  end
```

mRNA is produced at rate `k_y` from every active state and degrades at rate
`gamma` per molecule. Under the OR gate the active set is `{10, 01, 11}` for HD
and `{S, N}` for M — so the two topologies differ only in whether the factors
can be bound at the same time, which is exactly the question the counts are
asked to settle.

## Goals

1. **Predict, don't just simulate.** Give each architecture's stationary count
   distribution exactly — closed-form occupancies, burst frequency and size,
   and the Fano factor where they exist; finite state projection where they do
   not. Simulation is kept as a check on the algebra, not as the source of it.
2. **Decide between the architectures from data.** Fit both to allele-resolved
   counts for *sox2*, *nanog*, *rex1* and *esrrb* by sequential Monte Carlo —
   or fit them jointly, with the topology as a parameter, so the posterior
   reports P(heterodimer) directly. Validate the machinery against synthetic
   data generated from known parameters.
3. **Be honest about identifiability.** Stationary counts fix rates only in
   units of the mRNA degradation rate `gamma`, which is therefore held at 1 and
   never reported as inferred.
4. **Reproduce every figure from source**, with one script per figure and a
   single shared style module.

## Results at a glance

**Burst frequency per mRNA lifetime is the shape parameter.** `f/gamma` — how
many bursts arrive before a transcript decays — sets the shape of the count
distribution: large `f/gamma` gives a near-Poisson peak, small `f/gamma` gives
the long-tailed, over-dispersed distribution seen in the data. Everything else
about the promoter enters only through the burst size `b`.

![f over gamma as the negative-binomial shape parameter](figures/output/fig4_f_over_gamma.svg)

```bash
python figures/fig04_f_over_gamma.py
```

**How long the promoter stays on.** The mean first-passage time out of the
bound states, `tau_ON`, across the binding and unbinding rates of both factors.
The upper row varies one rate at a time; the lower row is the full plane, with
the parameters implied by single-molecule tracking marked.

![Mean first-passage time out of the bound states](figures/output/mfpt_alpha_beta.png)

The MFPT calculation lives in `figures/fig_1_mfpt_heatmaps.py`, which writes
`figures/output/fig_1_mfpt.svg`; the panel above is the tracked copy of the
same quantity.

## Install

```bash
git clone https://github.com/pacifistsilver/Srinjan_Modelling.git
```

```bash
cd Srinjan_Modelling && pip install -e .
```

Installing the package is what makes every script runnable from any working
directory. The core install covers the analytical results, the SSA, the CME
solver and the figures. The Bayesian inference stack is a heavier extra:

```bash
pip install -e ".[inference]"
```

```bash
pip install -e ".[all]"
```

### Docker

```bash
docker build -t stochtf .
```

```bash
docker run --rm -it -v $(pwd):/app stochtf python figures/fig04_f_over_gamma.py
```

## Quick start

Four commands, one per layer of the project:

```bash
python -c "from stochtf.analytical import heterodimer as hd; print(hd.fano(0.5, 0.05, 0.3, 0.2, k_y=30.0, gamma=1.0))"
```

```bash
python scripts/run_ssa.py --model heterodimer
```

```bash
python scripts/run_fsp.py heterodimer
```

```bash
python scripts/run_inference.py --gene esrrb --model heterodimer
```

The first prints the exact Fano factor with no simulation at all; the second
simulates one trajectory and writes a time-course plot to `results/`; the third
solves the master equation by finite state projection; the fourth fits the
model to real counts and writes a trace to
`results/esrrb_heterodimer.nc`.

## Usage

### Analytical results

Every quantity below is exact — no truncation, no Monte Carlo. Rates are in
units of `gamma`.

```python
from stochtf.analytical import heterodimer as hd

p = dict(a_s=0.5, b_s=0.05, a_n=0.3, b_n=0.2)   # SOX2 on/off, NANOG on/off

hd.pbound(**p)                        # P(at least one site occupied)
hd.burst_frequency(**p)               # f = 1 / (<T_on> + <T_off>)
hd.burst_size(**p, k_y=30.0)          # b = k_y <T_on>
hd.mean_y(**p, k_y=30.0, gamma=1.0)   # <y> = b f / gamma
hd.fano(**p, k_y=30.0, gamma=1.0)     # exact Fano factor, N-independent
```

The same four-state promoter under all three emission rules, mean and Fano side
by side:

```python
from stochtf.analytical import gates

for gate in ("OR", "AND", "ADD"):
    mean, fano = gates.analytic(0.5, 0.05, 0.3, 0.2, k=30.0, g=1.0, gate=gate)
    print(f"{gate:>3}  mean {mean:7.2f}  Fano {fano:6.2f}")
```

`gates.fsp(...)` returns the full stationary distribution for any of the three
gates and `gates.moments(P)` reduces it to `(mean, variance, Fano)` — the
independent check that the closed forms above are right.

Two results generalise past the two-site promoter: `analytical.dimer` gives the
exact Fano factor for *any* autonomous promoter-and-pool driver, bimolecular
dimerisation included, and `analytical.bursts` gives burst frequency and size
for any driver chain, verifying `<y> = b f / gamma` in each case.

### Simulate a trajectory

```bash
python scripts/run_ssa.py --model heterodimer
```

```bash
python scripts/run_ssa.py --model monomer --t-max 100 --out-dir /tmp
```

`--model` takes `monomer`, `homodimer` or `heterodimer`. Each writes a
time-course plot and prints the ON/OFF statistics extracted from the run.

### Solve the master equation

```bash
python scripts/run_fsp.py burr08 --expander support
```

```bash
python scripts/run_fsp.py heterodimer
```

The positional argument selects the network (`burr08`, `heterodimer`,
`homodimer`); `--expander` chooses between the `simple` and `support`
state-space expanders.

### Fit to data

```bash
python scripts/run_inference.py --gene esrrb --model heterodimer
```

```bash
python scripts/run_inference.py --gene sox2 --model joint
```

`--gene` takes `sox2`, `nanog`, `rex1`, `esrrb`, or any of the `synthetic_*`
datasets under `data/synthetic/` for validating the fit against parameters
known in advance. `--model` takes `monomer` (one contested site),
`heterodimer` (two independent sites), or `joint`.

`joint` is the interesting one: it makes the topology itself a parameter.
`model_index` is 0 for the monomer and 1 for the heterodimer, SMC moves
between them, and the posterior mean of the index *is* P(heterodimer | data).
Both topologies carry the same five rates under the same priors, so there is
no Occam asymmetry built into the comparison. Traces are written to
`results/<gene>_<model>.nc`.

The fit is scored against the exact stationary likelihood from
`stochtf.analytical.pgf` rather than by comparing a Gillespie run to the data
through a summary statistic — so the whole distribution is used and there is
no ABC tolerance to tune. Sampling is still SMC.

### Where the data can settle the topology

```bash
python scripts/run_inference.py --grid --stride 10 --draws 300
```

Runs the joint fit at every point of the synthetic grid in
`data/synthetic/heterodimer_grid.npz` and returns P(heterodimer) as a matrix
over `(alpha_s, alpha_n)`. Every dataset in that grid *was* generated by the
heterodimer, so the matrix maps where in parameter space the counts can
identify the topology at all: values near 1 mean the data settle it, values
near the 0.5 prior mean they cannot — which is the expected outcome wherever
the two promoters put nearly the same law on the counts.

A fit is not cheap (roughly 17 s at 200 draws on one chain), so the full 50×50
grid is about twelve hours. `--stride 10` gives a 5×5 map in minutes and
`--stride 5` a 10×10 in under an hour. Thin first, then decide whether the fine
grid is worth it.

### Figures

```bash
python figures/fig04_f_over_gamma.py
```

Each script writes an SVG into `figures/output/`. Two naming conventions are
present after the merge: `figNN_*.py` (the publication set mapped in
[figures/README.md](figures/README.md)) and `fig_N_*.py` (mean first-passage
times, burst kinetics, SSA time traces, distributions, model selection). The
numbering still needs confirming against the manuscript, and the two sets need
reconciling.

`figures/standalone_f_over_gamma.py` needs only numpy, scipy and matplotlib —
no install, no package import — and is the self-contained reference for the
`f/gamma` result.

### Notebooks

`notebooks/` holds the exploratory work: monomer noise (`01`), steady-state
binding (`02`), NANOG/SOX2 count histograms (`03`), CV² against the mean (`04`)
and posterior diagnostics for the fits (`05`).

## Data

Raw counts are allele-resolved single-cell RNA-seq from mouse ES cells, GEO
accession
[GSE132589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132589)
(Ochiai et al.). The two ~99 MB tables are downloaded rather than committed:

```bash
python scripts/download_data.py
```

Per-gene count vectors (`data/processed/*.npy`) and the synthetic validation
sets (`data/synthetic/*.npy`) are tracked, so the fits and figures run without
the download. See [data/README.md](data/README.md) for the full pipeline and a
known gap in it: `scripts/prepare_alleles.py` stops after building the
transcript-to-symbol map, so the processed arrays cannot currently be
regenerated end to end from the raw tables.

## Tests

```bash
pytest
```

```bash
pytest -m "not slow"
```

89 tests pass, with 4 xfail documenting known numerical discrepancies in the
analytical module.

`test_analytical.py` checks the identities the analytical module must satisfy:
`f = 1/(T_on + T_off)`, `<y> = b f / gamma`, the closed-form Fano factor
against exact FSP, and the single-site limits. `test_gillespie.py` checks the
simulator against them in turn. `test_pgf.py` and `test_stationary.py` check
the three routes to the stationary distribution against each other — closed
form, FFT inversion and direct ODE integration agree to machine precision.
`test_likelihood.py` and `test_identifiability.py` cover the exact likelihood
and the Fisher-information analysis; `test_models.py` needs the `inference`
extra installed.

## Layout

```text
├── src/stochtf/          installed package
│   ├── analytical/       closed-form + FSP results for two-site promoters
│   │                     (pgf, gates, heterodimer, monomer, dimer, bursts, mfpt)
│   ├── ssa/              Gillespie simulator, parameters, promoter models
│   ├── cme/              chemical master equation solver (FSP)
│   ├── inference/        exact likelihood, SMC models, identifiability
│   ├── plotting.py       shared figure style, palette, output paths
│   └── paths.py          project directory resolution
├── figures/              one script per paper figure -> figures/output/
├── scripts/              entry points (download, prepare, run, fit)
├── notebooks/            exploratory and diagnostic notebooks
├── tests/                pytest suite
└── data/                 processed and synthetic arrays (raw is downloaded)
```

## Citing

See [CITATION.cff](CITATION.cff). Please cite both the paper and the software.

## License

MIT — see [LICENSE](LICENSE).
