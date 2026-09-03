# PBS jobs for the complete certificate

For the complete paper calculation, use
`run_complete_proof_smoke.pbs`, then
`run_complete_proof_probe.pbs`, and finally
`run_complete_proof_production.pbs`.  These scripts are pinned to the
fixed logical staging directory
`/home/rendo/codex-jobs/Phanuel/20260903-104720-complete-proof/`; submit
the production job only after the smoke and probe jobs succeed.  Because `/home`
on the 2026-09-03 run was full, that logical job directory contains links
to the physical checkout and certificate workspace at
`/dev/shm/rendo-phanuel-20260903-104720-complete-proof/`, including the
checkout, certificate output, logs, MATLAB preferences, temporary files,
and local-pool job storage.  This node-local path is volatile and the final
artifacts must be retrieved immediately after validation.  All three jobs
must be submitted from the physical `state/` directory (the logical symlink
to it is also acceptable), because each wrapper verifies
`readlink -f "$PBS_O_WORKDIR"` before doing work:

```sh
cd /home/rendo/codex-jobs/Phanuel/20260903-104720-complete-proof/state
qsub ../source/hpc/run_complete_proof_smoke.pbs
# After the smoke completion marker and logs have been checked:
qsub ../source/hpc/run_complete_proof_probe.pbs
# After the probe summary and resource headroom have been checked:
qsub ../source/hpc/run_complete_proof_production.pbs
```

The wrappers take one shared kernel `flock`, so accidental concurrent runs
against the same checkout fail immediately.  The first wrapper atomically
records the staged Git commit; every later wrapper and production restart
requires the same `HEAD`.  Smoke and probe additionally require a clean
checkout.  The complete smoke requests two CPU slots because its final
runtime-binding test checks two independent MATLAB workers.  The Omega-mid
publication may intentionally modify only its
canonical result files during production, so a production restart checks
the pinned commit and leaves the driver's hash-bound publication journal to
validate or recover that state.

INTLAB stores proof-critical `INTLAB_CONST` data in a shared platform MAT
file and rewrites that file at every `startintlab` call.  All client and
worker initializations therefore take an interprocess Java NIO lock beside
the INTLAB root.  Runtime provenance hashes the complete INTLAB tree and
the complete platform MAT payload; it normalizes only the 24-character
`Created on` timestamp in the strictly validated MAT-v5 header.  The exact
relative filename, byte count, normalization policy, and normalized MAT
SHA-256 are included in the runtime record.  Platform text, MAT control
bytes, serialized constants, and every other INTLAB file remain bound.

The production script refuses to start unless at least 40 GiB of tmpfs
scratch and 600,000 free tmpfs inodes remain.  Because tmpfs consumes the
same physical memory, the memory and scratch budgets are not independent.
The two-worker probe is capped at 25 minutes and evaluates four
production-parameter Omega-mid
cells chosen to exercise direct CR, terminal LG, and adaptive subdivision.
It records GNU `time -v` (including maximum resident set size), tmpfs and
memory availability before and after the run, and checkpoint/output byte
sizes in `state/probe-summary-<PBS_JOBID>.txt`.  Its in-process assertions
check the exact cell IDs, interval mode, worker/chunk/checkpoint counts,
subdivision and FEM limits, mesh floors, runtime recheck, and zero errors.
The final `state/probe-complete.sha256` marker is written only after MATLAB
exits successfully and hashes the attempt record, summary, timing, logs,
resource samples, and every probe output.  Production recomputes all of
those hashes, validates the marker path set, and rejects a marker from a
failed later attempt.  This small sample is an
empirical sizing check, not a proof and not a guarantee of full-run peak
usage.  Production verifies the probe summary hash and source commit, scales
its measured aggregate two-worker RSS linearly to 20 workers, adds 25%, then
adds 25 GiB of anticipated tmpfs growth and a 16 GiB system reserve.  The
result (never less than 128 GiB `MemAvailable`) is the combined start gate;
an RSS estimate above 150 GiB is rejected as incompatible with the guarded
160 GiB PBS allocation.  The computed and observed values are preserved in
`state/production-sizing.txt`.  This makes the memory gate evidence-based
rather than fixing it above the memory actually available at submission,
while the independent 40 GiB free-space gate retains the storage-audit
safety margin.  The script sets
`VER10_PARALLEL_JOB_STORAGE`, so every local-pool worker writes its control
files to the job scratch rather than the full home filesystem.

The production allocation is one node, 20 CPU cores, 160 GiB RAM, and a
12-hour wall-time cap.  It creates the CR--LG spectral atlas and residual
certificate, the replayable thin-sector certificate, the archived affine
atlas validated from its separately clean immutable source commit, the
hash-bound global upper manifest, and finally the joint Omega-mid certificate.  The
middle-region raw MAT evidence is addressed by repository, release tag,
asset name, and SHA-256 rather than being committed to Git.
Once each merged certificate has passed its complete-certificate gate,
the production script removes only its exact expected per-task or per-chunk
checkpoint MAT set.  Cleanup requires the complete phase hash inventory,
the exact file count and canonical filename-set digest, no extra entry of
any type, and successful directory removal; otherwise it retains the
directory.  The final MAT, CSV, manifests, configs, and immutable
publication generation remain; the removed checkpoints are reproducible
from the pinned commit.

After a scheduler or MATLAB failure, resubmit the same production script
from `state/` without changing `source/`.  A valid relative-path SHA-256
inventory skips a completed upper or Omega-mid phase; an invalid inventory
forces validation/recomputation.  Generated hidden staging paths are never
silently incorporated into a completion marker: the wrapper stops and asks
the operator to inspect them before retrying.  Do not remove the volatile
scratch tree until its final marker, artifacts, and logs have been copied
off-node and rehashed.

The focused scripts below remain useful when only the thin Rayleigh
component needs to be regenerated.

The two scripts in this directory assume a fresh PBS job directory with a
clean checkout in `source/`.  Submit from the job directory, not from the
checkout:

```sh
qsub source/hpc/run_muhat1_rayleigh_smoke.pbs
qsub source/hpc/run_muhat1_rayleigh_full.pbs
```

Both jobs use queue `normal`, one CPU, 4 GiB of memory, and single-threaded
BLAS.  The smoke job runs the focused interval replay/tamper test and has a
30-minute wall limit.  The production job has a one-hour wall limit and
executes exactly

```sh
matlab -batch "run('scripts_run/run_muhat1_rayleigh_hpc.m')"
```

By default, the production artifacts are written outside the checkout to
`$PBS_O_WORKDIR/certificate-output/muhat1_rayleigh_t038/`:

- `manifest.json`;
- `frozen_coefficients.csv` (lossless `q` and final monomial intervals);
- `cells.csv` (all Rayleigh cell enclosures);
- `down_cells.csv` (all Cartesian `(s,r)` cell enclosures);
- `checksums.sha256`.

The default INTLAB location is the verified Liu Lab copy
`/home/rendo/tmp/Intlab_V12_nqc_01`.  Override `INTLAB_ROOT`,
`THIN_SOURCE_DIR`, or `VER10_CERTIFICATE_ROOT` in the `qsub` environment if
the staged layout differs.  The full driver rejects a dirty source tree,
an in-tree output path, noncanonical parameters, existing artifacts, hash
mismatches, and any nonnegative replayed upper endpoint.
