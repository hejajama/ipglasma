# Merging `optimized_bjorn` into `merge_performance_improvements`

Notes on resolving the merge of `optimized_bjorn` (Bjoern Schenke's
performance rewrite: SoA `Lattice` matrix storage, SU(3) specialization,
`Instrumentation` profiling, bulk-RNG `GaussBulk`, ...) into this branch,
which carries JIMWLK small-x evolution (`jimwlk.cpp`, `mode`-based `Cell`
observables) that `optimized_bjorn` never had.

Conflicted files: `input`, `src/Cell.cpp`, `src/Cell.h`, `src/Init.cpp`,
`src/Init.h`, `src/Lattice.cpp`, `src/Lattice.h`, `src/Random.cpp`,
`src/Random.h`, `src/main.cpp`.

## Resolution direction

Adopted `optimized_bjorn` wholesale for `Cell`, `Lattice`, `Random`: matrices
moved off `Cell` onto flat `std::vector<Matrix>` fields on `Lattice`
(`U`, `U2`, `Ux`, `Uy`, `Ux1`, `Uy1`, `Ux2`, `Uy2` — logical aliases `U/E1`,
`U2/E2`, `Ux1/g`, `Uy1/Uplaq`, `Ux2/pi`, `Uy2/phi`), `Cell` now holds only
scalar per-site observables (`Tmunu`/`pimunu`/`umu` components as individual
doubles, no `mode_` flag).

That removed the `Cell::getU()/setU()`-style API that `jimwlk.cpp` and
`Lattice::WriteWilsonLines`/`WriteSU3Matricies` (JIMWLK-only, absent from
`optimized_bjorn`) depended on, and `optimized_bjorn`'s rewritten
`Init::init()` fused impact-parameter sampling directly into
color-charge-density construction (see the "positions are shifted here, not
later" comment in the original) — which conflicts with JIMWLK's requirement
that each nucleus's Wilson lines be evolved *before* any impact parameter is
applied.

## What was ported / rewritten (not just re-pointed at new APIs)

- **`jimwlk.cpp`**: `Cell::getU()/setU()`/`getU2()/setU2()` → `Lattice::U[i]`
  / `Lattice::U2[i]`. Mechanical, no behavior change.
- **`Lattice::WriteWilsonLines` / `WriteSU3Matricies`**: re-added (didn't
  exist on `optimized_bjorn`), ported to the SoA fields.
- **`Init::init()`**: split into phases so JIMWLK can run in between:
  - `init()` now builds each nucleus's color-charge density and Wilson
    lines in its own centered frame — `b` is forced to `0` (and
    `useFixedNpart`'s resampling disabled) for this call only, then
    restored.
  - Four new methods, ported from this branch's original (pre-merge)
    `Init.cpp`, apply the real collision geometry afterward:
    `sampleImpactParameter`, `computeCollisionGeometryQuantities`,
    `shiftFieldsWithImpactParameter` (ported to the new `BufferLattice` SoA
    API), `initializeForwardLightCone` (extracted verbatim from the tail of
    `optimized_bjorn`'s monolithic `init()` — that part was already
    SoA-ported and already improved over this branch's version, e.g. a
    completed "plus ax,ay" branch that was previously commented out, and
    per-cell retry seeds for determinism).
  - `main.cpp`'s event loop: build nuclei → JIMWLK evolution (if
    `useJIMWLK`) → `sampleImpactParameter`/`computeCollisionGeometryQuantities`
    → `shiftFieldsWithImpactParameter` → `initializeForwardLightCone` →
    `evolution.run()`.
- **`Parameters.h`**: added `writeInitialWilsonLines` (get/set), kept as a
  separate flag from JIMWLK's `writeWilsonLines` — they gate genuinely
  different dumps (initial condition vs. periodic JIMWLK snapshots). Kept
  the existing `rapidityA_`/`rapidityB_` split (not `optimized_bjorn`'s
  single `rapidity`) since `WriteWilsonLines`'s binary header needs
  per-nucleus rapidity; `main.cpp` parses `RapidityA`/`RapidityB` from
  `input` accordingly.
- **`Random.h`/`.cpp`**: kept both `optimized_bjorn`'s `GaussBulk` (bulk RNG)
  and this branch's `sampleGammaInc`/`setGammaIncCDF`/`ranGen_` — see caveat
  below, the latter is currently unused.
- **`input`**: kept `optimized_bjorn`'s values for the two conflicting
  blocks (`bmax`, `computeGluonMultiplicity`, `writeOutputs`, etc.), added
  back the JIMWLK parameter block and `writeWilsonLines`.

## Verification done

- Clean build (`cmake` + `make`, AppleClang, one pre-existing unrelated
  warning).
- End-to-end smoke test (small lattice, `useJIMWLK 1`): initialization →
  JIMWLK evolution on both nuclei → impact-parameter sampling
  (`Npart`/`Ncoll`/`Qs` all sane) → field shift → forward lightcone → full
  CYM evolution → gauge fixing converged (residual ~1e-9) → multiplicity
  written. No crashes. Gauss-law violation after evolution: ~7.7e-11.

## Known gaps / things that need a physics review, not just a compile check

1. **The non-JIMWLK path changed too.** Every event, JIMWLK or not, now goes
   through "build at `b=0`, then re-index the lattice by the sampled `b`"
   (`shiftFieldsWithImpactParameter`) instead of Bjoern's single-pass
   "bake `b` into the color-charge loop directly." These are two different
   discretizations of the same physical shift. The smoke test looked sane
   (tiny Gauss-law violation, reasonable `Npart`/`Qs`/`dN/dy`), but there
   has been no direct numerical comparison against an unmodified
   `optimized_bjorn` build for a non-JIMWLK event.
2. **`useFixedNpart` + `useJIMWLK` together is unsupported.** Bjoern's
   Npart-matching retry loop needs the real `b` baked into
   color-charge-density construction, which conflicts with building nuclei
   before `b` is known. Noted in a comment in `Init::init()`.
3. **`useFatTails`/`tDistNu` (fat-tailed proton position sampling) is gone.**
   `optimized_bjorn` dropped this feature entirely, independent of JIMWLK.
   `Random::sampleGammaInc()`/`setGammaIncCDF()` were kept (harmless,
   currently unused) rather than re-wiring them into Bjoern's rewritten
   `setColorChargeDensity`, since that was out of scope for the JIMWLK
   surgery. `input` currently has `useFatTails 0`, so this is silent unless
   turned on.

## Suggested follow-up

- Run a non-JIMWLK event through both an unmodified `optimized_bjorn` build
  and this merged branch with the same seed; compare `Npart`, `Ncoll`,
  `Qs`, and `dN/dy` to bound the effect of the shift-based geometry change.
- Decide whether `useFatTails` is worth re-porting into the new
  `setColorChargeDensity`, or should stay dropped.
- Decide whether `useFixedNpart` + `useJIMWLK` needs real support, or
  whether an explicit runtime check/error should reject that combination
  instead of silently ignoring `useFixedNpart`.
