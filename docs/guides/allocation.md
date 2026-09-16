# Allocation and memory in prakash

> **prakash never frees.** Every byte a call allocates is held for the life of the
> process. This page says what each entry point costs so a consumer can size its
> own usage, and what the two reclaim options are.

## The allocator

prakash allocates through the Cyrius stdlib's bump allocator (`lib/alloc.cyr`):
a pointer walks forward through mmap'd chunks, and `free` does not exist. There
is no garbage collector and no reference counting.

⚠ **`alloc()` does not take a lock in a single-threaded process.**
`_alloc_lock_acquire` returns immediately while `_threads_active == 0`, and
prakash never spawns a thread. That flag is process-wide, so a *consumer* that
spawns one arms the lock for every subsequent allocation.

## Per-call cost of the public entry points

Measured on x86-64 with `alloc_used()` deltas, steady state (caches warm) at
2.3.8. "Retained" is the part you keep — the returned object. "Scratch" was
working memory; since 2.3.8 the 2D diffraction path reuses it.

| Entry point | Bytes per call | Note |
|---|---:|---|
| `trace_surface` (hit) | 128 | 80 on TIR; **0** on a miss or aperture reject |
| `trace_sequential` (2 surfaces) | 408 | 2 × `trace_surface` + the vec |
| `spot_diagram` (3 × 6 fan) | 11,304 | 38 × `trace_surface` |
| `rgb_to_json` | 136 | |
| `medium_to_json` | 136 | |
| `sellmeier_to_json` | 336 | |
| `prescription_to_json` (6 surfaces) | 1,480 | grows with surface count |
| `spd_to_json` (81 samples) | 2,576 | grows with `Spd_len` |
| `diffraction_pattern_2d(64×64)` | 32,792 | was 99,352 before 2.3.8 |
| `diffraction_pattern_circular(32)` | 8,216 | was 33,304 |
| `psf_from_wavefront(32×32)` | 8,216 | was 25,112 |
| `pbr_integrate_brdf_lut` | 0 | since 2.3.7 |
| `multilayer_rt` | 48 | 80 before 2.3.7 |

**Formulas for the pattern functions**, exact at every size measured:

```
diffraction_pattern_2d(n, n)      = 8n² + 24     retained
diffraction_pattern_circular(n)   = 8n² + 24     retained
```

⚠ **The returned `Pattern2D` is the floor and it is never reclaimed.** A
256-sample PSF sweep at 64×64 retains about 8.4 MB no matter what the scratch
does. Size for the retained figure, not the per-call one.

## Reclaiming memory

There are two options and **one of them is a trap**.

### 1. Process lifetime

The intended model. Load, compute, exit. No action needed.

### 2. `alloc_reset()` — and you MUST follow it with `prakash_reset_caches()`

⛔ **Calling the stdlib's `alloc_reset()` without `prakash_reset_caches()`
silently corrupts every colour value prakash produces.** prakash memoises 21
tables — the CIE observers, the standard illuminants, the CRI test-colour
samples, the photopic and scotopic V(λ) curves — as pointers into the heap that
`alloc_reset()` scrubs. The pointers survive; the data does not.

Measured before the fix: `cie_cmf_at(555)` returned ȳ = **1.0** before a reset
and **0.0** after. 555 nm is the photopic peak, which is 1.0 by definition, so
`spd_to_xyz`, `luminous_flux` and the whole CRI computation silently produced
garbage. No error was reported and nothing crashed.

```cyrius
alloc_reset();
var _r = prakash_reset_caches();   # ⛔ not optional
```

`prakash_reset_caches()` is a *forget*, not a free: it drops the memo pointers so
the next call rebuilds each table from the compiled-in constants. Rebuilding all
21 costs about 24 KB.

## Thread safety

prakash makes **no thread-safety guarantee**, and two things narrow it further:

- The **2D diffraction path** (`diffraction_pattern_2d`,
  `diffraction_pattern_circular`, `psf_from_wavefront`,
  `psf_diffraction_limited`) uses module-scope scratch caches since 2.3.8.
  Concurrent calls into it will corrupt each other. Serialise them.
- The 21 memoised tables are built lazily. Two threads racing the first call
  build identical bytes, so that race is benign — but `prakash_reset_caches()`
  during a concurrent call is not.

Everything else allocates per call and holds no shared state.

## If you are running out of memory

1. Are you retaining `Pattern2D` results you no longer need? That is the floor,
   and no prakash change can lower it.
2. Can you restructure as compute-then-exit? That is the model prakash is for.
3. If neither, `alloc_reset()` + `prakash_reset_caches()` between epochs is the
   only reclaim, and it invalidates every handle you still hold.
