"""Parameter inference against allele-resolved count data.

Modules:
    models: SMC over the two competing promoter topologies -- the heterodimer
      (two independent sites) and the monomer (one contested site) -- or over
      both at once, with the topology itself as a parameter.
    likelihood: The exact stationary log-likelihood, built on the distribution
      from ``stochtf.analytical.pgf``. Covers the independent-site gates only.
    identifiability: What stationary counts can determine at all, from the
      Fisher information. Worth reading before interpreting any marginal.
"""
