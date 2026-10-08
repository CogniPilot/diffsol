# Rumoca 0.13.2 maintenance branch

This branch preserves upstream commit
`7edd02af717287394bac346d4ed74ec141ae5a86` and Diffsol's published 0.13.2
version and API. Original authorship and MIT licensing are unchanged.

The changes in `diffsol/src/ode_solver/bdf.rs` and
`diffsol/src/nonlinear_solver/line_search.rs` are copied byte-for-byte from
Rumoca's former vendored dependency. They publish the corrected accepted BDF
endpoint, enforce the configured minimum before an attempt while respecting
hard stops, refresh the Newton coefficient after a shortened hard stop,
consume accepted Newton correction history and apply pending corrections, and
scope nonlinear failure limits to the current step while retaining cumulative
statistics. Tolerances, Armijo conditions, and failure limits are unchanged.

`diffsol/tests/bdf_numerical_regressions.rs`, `convergence_tests.rs`, and
`failure_budget_tests.rs` reproduce these numerical contracts independently of
Modelica and Rumoca's event host. They originated in
`CogniPilot/rumoca` at `2aaed750242f96d77f39e319f16b24bee82794c1`.

Run the focused suite with:

```sh
cargo test -p diffsol --test bdf_numerical_regressions \
  --test convergence_tests --test failure_budget_tests
```

This is a maintenance fork, not an upstream release. It does not claim that
upstream main has these corrections, nor that it has completed Rumoca's full
MSL/OMC qualification. Rumoca records the complete patch provenance and broader
validation requirements in its simulation runtime documentation.
