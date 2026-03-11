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

### `phydrus/read.py` — Numeric Index & Data Preservation Fix

A critical bug was found in `read_obs_node`:
1.  **Unconditional Row Dropping:** The parser was unconditionally dropping the first row of data, assuming it always contained units (e.g., `[T]`). In simulations with low print frequency (like Example 1 in the notebook), this often **deleted the only data point**, resulting in empty DataFrames.
2.  **String Indices:** When units rows *were* present, the `time` index was parsed as strings. Plotting strings against numeric axes (e.g., `xlim(0, 730)`) resulted in **empty plots**.

**Fix:**
- Updated `read_obs_node` and `_read_file` to **robustly detect units rows** by checking for brackets `[` or `]`.
- Enforced **numeric conversion** for both the DataFrame columns and the index (`pd.to_numeric`).
- Verified that **single-time-point data** (e.g., at time 730) is now preserved and plottable.

---

## Final Verification Results

- **12 passed** in the core test suite.
- **`example_1.py`**: Successfully generates Bottom Flux and Water Content plots.
- **`example_2.py`**: Successfully completes long-running simulation loop.

Verified with diagnostic plotting script:
![Water Content Plot](/Users/veethahavya/.gemini/antigravity/brain/6bf7f9e9-b7b7-41ce-8f7b-d8cd712d89e8/test_water_content.png)
![Bottom Flux Plot](/Users/veethahavya/.gemini/antigravity/brain/6bf7f9e9-b7b7-41ce-8f7b-d8cd712d89e8/test_water_flow.png)
