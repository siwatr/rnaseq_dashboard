# Object states, edit history, and the reset/undo controls

Developer reference for how the app tracks the dataset and what each "reset"/"undo"
affordance actually does. The canonical store is the app-state object in
[R/state.R](../R/state.R); see also CLAUDE.md "State model".

## The three object states

- **`original`** — the dataset exactly as loaded. **Never mutated.** It is the target
  for every "reset to original".
- **`working`** — the current `dds` after committed edits (metadata, annotation, assays,
  sample/feature removal). This is what every page reads and plots.
- **`draft`** — a *local, uncommitted* buffer that lives **only inside the metadata
  editor** (`mod_meta_editor`). Cell edits / add-column / etc. accumulate here and only
  reach `working` when the user clicks **Save** (one `state_mutate`).

Every committed edit goes through `state_mutate()`, which: pushes the pre-edit `working`
onto **`undo_stack`** (LIFO, hard-capped at **depth 5** — `.undo_depth`), applies the edit,
clears the `derived` cache, bumps **`data_version`**, and appends to the `history` log.

```
  original --load--> working --edit1--> edit2 --> edit3(filter) --> [ current ]
  (immutable;                                                            |
   the reset target)                                                     |
                                                                         |
  per-edit snapshot pushed onto undo_stack (LIFO, max depth 5) ----------+

  GLOBAL  (status bar, always visible)
    [Undo ]  step back the LAST committed edit only  (any kind; whole snapshot)
    [Reset]  working <- original                     (drops ALL edits)

  LOCAL  (inside a tab; scoped - unrelated edits untouched)
    metadata "Reset to original"   revert ONLY this slot (colData | rowData)
                                   to original values, current rows
    QC "Reset Sample/Feature       re-add ONLY the removed samples | features
        Removal"                   from original; keep edits on kept items
                                   and the other slot

  draft (metadata editor only): a local, UNCOMMITTED buffer.
    Save               -> commit draft via state_mutate (enters the line above)
    Reset to last save -> draft <- working   (throw away uncommitted edits)
    "Unsaved changes" badge shows when draft differs from working.
```

## The controls, by scope

| Control | Where | What it reverts | Scoped? | Function |
|---|---|---|---|---|
| **Undo last edit** | status bar (global) | the **single most recent** committed edit, as a whole-object snapshot (LIFO, depth 5) — could be a metadata edit made *after* a filter | **No** (order-sensitive) | `state_undo()` |
| **Reset to original** | status bar (global) | `working <- original`; **all** edits; clears undo history | No (full) | `state_reset()` |
| **Reset Sample/Feature Removal** | QC Filtering pill (local) | re-adds the removed samples *or* features from `original`; keeps edits on kept items **and the other dimension/slot** | **Yes** (one dimension) | `restore_samples()` / `restore_features()` |
| **Reset to original** | metadata editor (local) | that one slot (`colData` *or* `rowData`) to original values for the current rows; keeps filtering, assays, the other slot | **Yes** (one slot) | `reset_metadata_slot()` |
| **Reset to last save** | metadata editor (local) | the uncommitted `draft` back to `working` (pre-commit only) | n/a (draft) | `draft(working)` |

Key consequence: the **local** resets are the surgical ones — e.g. "Reset Sample Removal"
re-adds removed samples but does **not** touch rowData edits or feature-filtering. The
**global Undo** is *not* slot-scoped: it reverts whatever the last committed action was,
because it restores a whole snapshot. There is no per-tab "undo one step" — locally you get
targeted **resets** only.

## Implementation notes

- All four reset/restore paths are themselves committed via `state_mutate`, so each is
  **undoable** (it lands on `undo_stack`) and appears in `history`.
- `restore_*` / `reset_metadata_slot` **merge metadata**: items that stayed keep their
  edited values; re-added items take `original`'s, aligned to the *current* column schema
  (columns added after load become `NA`). Counts come from `original` (they are never
  edited); normalized assays are recomputed (`refresh_assays`) and size factors
  re-estimated when they were set.
- The metadata "Reset to original" commits immediately (no Save) and resyncs the `draft`
  via the editor's `observeEvent(state$working, …)`.
- Memory: `undo_stack` holds at most 5 prior `dds` snapshots. R shares unchanged memory,
  but subsetting/edits allocate fresh assay matrices, so the worst case is ~5x the
  mutated-assay footprint — bounded, and fine for bulk; lower `.undo_depth` if a future
  very-large/single-cell path needs it.
