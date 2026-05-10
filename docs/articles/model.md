# Model structure

This page sketches the methodology behind `ineqx`. Three equations carry
the intuition; for the full derivations see [Rosche
(2026)](https://osf.io/preprints/socarxiv/f53kz/).

### 1. Descriptive within/between decomposition

Take the vector $`Y_t`$ of individual incomes at time $`t`$, and the
vector $`G_t`$ of group memberships ($`G_t = j`$ for $`j = 1,\dots,J`$
mutually exclusive groups). The variance of $`Y_t`$ splits cleanly into
a within- and a between-group part:

``` math
\begin{equation}
V_t \;=\; \underbrace{\sum_j \pi_{jt}\,\sigma_{jt}^2}_{\text{within } W_t}
       \;+\; \underbrace{\sum_j \pi_{jt}\,\big(\mu_{jt} - \bar\mu_t\big)^2}_{\text{between } B_t}
\tag{1}
\end{equation}
```

with $`\pi_{jt}`$ the share of group $`j`$, $`\mu_{jt}`$ its mean,
$`\sigma_{jt}^2`$ its variance, and
$`\bar\mu_t = \sum_j \pi_{jt}\mu_{jt}`$ the overall mean. Cross sections
collapse to a single $`t`$. With repeated cross-sections, the change
$`V_t - V_{t_0}`$ further decomposes into contributions from changing
means $`(\Delta_\mu)`$, changing dispersions $`(\Delta_\sigma)`$, and
changing group shares $`(\Delta_\pi)`$, with each component split into
its between- and within-group share. This is what
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md)
returns when no `treat` is supplied — it is essentially the Western &
Bloome (2009) decomposition with a Shapley-averaged path.

### 2. Causal cross-sectional decomposition

Let $`D \in \{0,1\}`$ be a binary treatment, $`\beta_j`$ its causal
effect on group $`j`$’s mean, and $`\lambda_j`$ its causal effect on
group $`j`$’s log-SD (so the treated SD is
$`\sigma_j(0)\cdot e^{\lambda_j}`$). The treatment effect on inequality,
$`\tau = V[Y(1)\mid D] - V[Y(0)\mid D]`$, splits into a between- and a
within-group part, **and each of those further splits into a
heterogeneity and a covariance term**:

``` math
\begin{equation}
\boxed{\;\tau_B \;=\; \underbrace{\mathrm{Var}_\pi(\beta)}_{\text{het}_B}
            \;+\; \underbrace{2\,\mathrm{Cov}_\pi(\mu(0),\,\beta)}_{\text{cov}_B}
\;}
\qquad
\boxed{\;\tau_W \;\approx\; \underbrace{2\,\bar\sigma^2(0)\,\bar\lambda}_{\text{het}_W}
            \;+\; \underbrace{2\,\mathrm{Cov}_\pi(\sigma^2(0),\,\lambda)}_{\text{cov}_W}\;}
\tag{2}
\end{equation}
```

The four sub-components have direct substantive readings:

- **$`\text{het}_B = \mathrm{Var}_\pi(\beta) \ge 0`$.** Heterogeneity in
  the treatment effect on the mean *always* widens between-group
  inequality. Equal $`\beta`$ across groups means $`\text{het}_B = 0`$.
- **$`\text{cov}_B = 2\,\mathrm{Cov}_\pi(\mu(0),\beta)`$.** Sorting of
  gains. If positive, larger gains accrue to already-advantaged groups
  (reinforces the hierarchy). If negative, treatment compresses the
  hierarchy.
- **$`\text{het}_W \propto \bar\lambda`$.** The average effect on
  within-group dispersion. Sign tracks $`\bar\lambda`$: treatment can
  stretch or compress within-group spread.
- **$`\text{cov}_W = 2\,\mathrm{Cov}_\pi(\sigma^2(0),\lambda)`$.**
  Sorting of variance effects. If positive, treatment amplifies
  dispersion most in already-disperse groups.

`plot(<causal>, type = "wibe", stats = c("het_b","cov_b","het_w","cov_w"))`
shows these four series directly.

The package supports both $`V`$ (variance) and $`CV^2 = \sigma^2/\mu^2`$
(squared coefficient of variation) as the underlying inequality measure.

### 3. Longitudinal three-channel decomposition

With repeated cross-sections, the *change* in the treatment effect on
inequality, $`\Delta\tau = \tau_t - \tau_{t_0}`$, decomposes into three
interpretable channels:

``` math
\begin{equation}
\Delta\tau \;=\; \underbrace{(\Delta_\beta + \Delta_\lambda)}_{\text{behavioral}}
\;+\; \underbrace{(\Delta_{\pi,B} + \Delta_{\pi,W})}_{\text{compositional}}
\;+\; \underbrace{(\Delta_\mu + \Delta_\sigma)}_{\text{pre-treatment}}
\tag{3}
\end{equation}
```

- **Behavioral** — changes in *how* treatment affects each group’s mean
  and variance ($`\Delta_\beta`$, $`\Delta_\lambda`$). The
  Kitagawa-Blinder-Oaxaca “coefficient effect.”
- **Compositional** — changes in *who gets treated* across groups
  ($`\Delta_{\pi,B}, \Delta_{\pi,W}`$). The KBO “endowment effect.”
- **Pre-treatment** — changes in the baseline distribution that
  determine how much a fixed treatment effect *translates into*
  inequality ($`\Delta_\mu, \Delta_\sigma`$). Has no KBO analogue — this
  is what the variance-decomposition framing adds.

Sequential decompositions of this kind are path-dependent: the size of
$`\Delta_\beta`$ depends on whether you switch $`\beta`$ before or after
$`\pi`$.
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md)
defaults to `order = "shapley"`, which averages across all $`3! = 6`$
orderings to give an order-invariant decomposition.
`plot(<causal>, type = "shapley")` displays the Shapley point estimate
and the range across orderings.

### 4. A note on the variance of logs

[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md)
accepts `ystat = "VL"` for the descriptive case (it log-transforms $`y`$
and runs the standard variance decomposition), but the paper recommends
**against** $`V_L`$ for substantive decomposition work. The reason:
because the log transformation is non-linear, the within- and
between-group components of $`V_L`$ can move in the *opposite direction*
from within- and between-group changes on the income scale. A treatment
that unambiguously raises within-group dispersion in dollars can
register as compressing it in $`V_L`$. $`V`$ and $`CV^2`$ do not have
this pathology and decompose into components that remain interpretable
on the original income scale.
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md) issues
a runtime warning whenever `ystat = "VL"` is used — see the FAQ for
details.
