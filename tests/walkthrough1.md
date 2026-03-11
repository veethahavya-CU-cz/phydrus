# Phydrus Python 3.14 Port — Walkthrough

## Root Cause: Pandas 3.0 MultiIndex Rendering Change

The `self.materials` DataFrame in `phydrus` uses a **MultiIndex** for columns, where the top-level group names are `'water'` and `'solute'`. In **Pandas 3.0**, calling `to_string(index=False)` on a MultiIndex DataFrame now renders the top-level group label as a **separate header row** — something older Pandas versions didn't do.

Hydrus Fortran's `SELECTOR.IN` parser expects:
```
  thr  ths  Alfa    n    Ks   l
0.078 0.43 0.036 1.56 24.96 0.5
```

But Pandas 3.0 was generating:
```
water                              ← EXTRA LINE — Fortran crashes here
  thr    ths   Alfa      n   Ks   l
  0.0 0.3382 0.0111 1.4737 13.0 0.5
```

This caused the Fortran error: `Error when reading from an input file Selector.in Water Flow Informations!`. Because `T_LEVEL.OUT` was never populated (Hydrus crashed silently), plotting then failed with `ValueError: Could not find start string 'rTop'`.

---

## Changes Made

### `phydrus/model.py` — `write_selector()`

**Water materials block (Block B):** Strip MultiIndex top-level before `to_string()`:
```python
mat_df = self.materials.fillna(0)
if isinstance(mat_df.columns, MultiIndex):
    mat_df = mat_df.copy()
    mat_df.columns = mat_df.columns.get_level_values(-1)
lines.append(mat_df.to_string(index=False, justify="right"))
```

**Solute materials block (Block F):** Filter to only `'solute'` columns (Fortran Block F expects ONLY `bulk.d, DisperL, frac, mobile_wc`):
```python
_full_mat = self.materials.fillna(0)
if isinstance(_full_mat.columns, MultiIndex) and "solute" in _full_mat.columns.get_level_values(0):
    _sol_mat = _full_mat["solute"].copy()
...
lines.append(_sol_mat.to_string(index=False, justify="right"))
```

**Per-solute data (`sol['data']`):** Same MultiIndex stripping applied to each solute's parameter DataFrame.

### `phydrus/read.py` — `_read_file()`

Made the start/end line search robust:
- Initialized `s = None`, `e = len(lines)` before the loop
- Only sets `s` once (first match), then looks for `e` only after `s` is found
- Raises a clear `ValueError` if start string is never found

### Example Scripts

Fixed all relative path issues so examples work regardless of working directory:
- **`example_1.py`**: `ws`, `exe`, ATM CSV reads, and CSV writes all use `os.path.join(_script_dir, ...)`
- **`example_2.py`**: Same fixes for `ws`, `exe`, and ATM CSV read

---

## Test Results

```
12 passed in 0.71s
```

All tests pass:
- `test_001.py::test_import`
- `test_read.py` — 9 output parsing tests (T_LEVEL, NOD_INF, OBS_NODE, RUN_INF, I_CHECK, PROFILE, A_LEVEL, BALANCE, SOLUTE)
- `test_write_input.py` — NaN-free SELECTOR.IN and PROFILE.DAT generation

## Example Runs

Both `example_1.py` and `example_2.py` now run end-to-end:
- `example_1.py`: Water content at obs nodes + Bottom Flux plot + profile update — ✅
- `example_2.py`: Multi-run solute transport simulation (250+ Hydrus runs) — ✅
