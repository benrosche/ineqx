# Installation

### Installation consists of two steps:

1.  Install devtools if necessary: `install.packages("devtools")`
2.  Install ineqx from GitHub:
    `devtools::install_github("benrosche/ineqx")`

The current release is **v0.4.0**. See
[NEWS.md](https://github.com/benrosche/ineqx/blob/master/NEWS.md) for
the changelog (v0.4.0 introduced the unified
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md) entry
point and is a breaking change relative to v0.3.x).

### Notes:

- `ineqx` requires `gamlss` for the causal estimation path. If GAMLSS is
  not yet installed, R will prompt for it on first use.
- In some cases, dependencies must be installed or updated before
  `ineqx` can be installed.
- If you get `Error: (converted from warning) ...`, you can set
  `Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")`
- For questions or feature requests, please open an issue at
  [github.com/benrosche/ineqx/issues](https://github.com/benrosche/ineqx/issues).
