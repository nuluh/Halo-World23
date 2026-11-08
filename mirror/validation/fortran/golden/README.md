# Golden Fortran oracle fixtures

This directory turns the raw deterministic oracle snapshots under `cases/3dframe/expected/` into a provenance-bearing fixture contract.

## Contract

`catalog.json` records, for every fixture:

- case identifier and physical/model inputs;
- legacy routine names and SHA-256 hashes of the untouched source files;
- reference compiler/version;
- driver, runner, and exact regeneration command;
- SHA-256 of the driver, runner, and frozen output CSV;
- a unit rule that covers every numeric key in the CSV;
- the tolerance policy used by downstream C++ parity tests.

The numeric values remain in the existing key/value CSV snapshots. They are not duplicated into JSON, so there is only one authoritative expected-value payload.

Run the metadata/integrity check without a Fortran compiler:

```bash
python3 validation/fortran/golden/verify.py
```

Regenerate all committed stage/progression fixtures from the preserved Fortran source and byte-compare them with the committed snapshots:

```bash
python3 validation/fortran/golden/verify.py --regenerate
```

The regeneration path uses the existing `run_*.sh` drivers and therefore honors `FC` and `VFIFE_ORACLE_WORKDIR`.

The exact GNU Fortran compiler and version used to freeze each fixture are stored in its `generation` record. The first stage-level fixtures were accepted with 14.2.0; the progressive frame/transient fixture added by #31 was frozen with 14.4.0 after the original stage fixtures reproduced under the same run.

## Tolerance policy

The stage-level oracle evidence produced by issues #28 and #29 is deterministic. The isolated central-difference equations and the untouched full `VFIFE()` routine measured `max_abs_diff = 0.0`.

For numeric unit class `u`:

```text
rtol = 256 * DBL_EPSILON
atol[u] = rtol * max(1, max_abs_expected[u])

pass iff
abs(actual - expected) <= atol[u] + rtol * abs(expected)
```

The factor 256 is an explicit floating-point roundoff budget. The absolute floor is tied to the largest magnitude actually present in that unit class, so a near-zero expected value is not judged by relative error alone.

This is intentionally a **stage-level baseline**, not permission to loosen later integration tests. If #31 or #32 measures reproducible cross-implementation drift above this baseline, update the policy with the measured distribution and cite the new evidence. Do not silently increase tolerances to make a test pass.

`central_step_oracle.csv:max_abs_diff` is dimensionally mixed because it takes one maximum across displacement, velocity, acceleration, and increment. It is retained as a diagnostic-only legacy metric; typed per-DOF errors are authoritative.

## Evidence boundary

The tolerance study is explicitly split by stage in `catalog.json`.

Measured now:

- geometry/local basis;
- lumped mass and rotational inertia;
- element forces;
- global force assembly;
- one explicit integration step;
- short transient history — measured by #31 against the untouched four-step VFIFE loop (72 values, maximum absolute error 5.087412011646086e-13);
- longer end-to-end regressions — measured by #32 with exact GNU Fortran 14.4.0 output signatures for the complete 30,000-step EX4 case and a deterministic 1,000-step prefix of L10E40.

The end-to-end regression criterion is exact signature identity, not a widened floating-point tolerance. Modern C++ drift is checked independently against the stage/progression fixtures. See `validation/fortran/regression/README.md` for the bounded L10E40 transform and the path to the original 2,000,000-step run.
