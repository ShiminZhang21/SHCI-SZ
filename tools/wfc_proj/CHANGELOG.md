# wfc_proj fix notes

## 2026-08-31 — `Compute_Lambda`: non-contiguous `A_S` / positional-slicing bugs

`Compute_Lambda`'s `G_S` membership test used to define `A_S` by a fixed
`AS_first`/`AS_last` block boundary (`index_map.front()`/`.back()`) and slice
`up_orbs`/`dn_orbs` by **array position** to check that the "core" was full
and the "virtual" was empty. This assumed `A_S` was a contiguous range of
`A_L` orbital indices, and broke silently in three ways found while
reviewing and testing it:

1. **Non-contiguous `A_S`.** A gap orbital — active in `A_L` but excluded
   from `A_S`, sitting between `AS_first` and `AS_last` — was invisible to
   the check: an electron occupying it was neither counted against `A_S` nor
   flagged as a leak, so a determinant with a genuinely misplaced electron
   could be wrongly admitted to `G_S`.
2. **`AS_last` vs. electron count.** The right-hand (virtual) guard compared
   an orbital *index* (`AS_last`) against an electron *count*
   (`up_orbs.size()`). Whenever `AS_last` exceeded that count — routine any
   time `A_L` has more active orbitals than electrons in that spin channel —
   the guard silently fell back to an empty slice, i.e. "no leakage
   detected," regardless of the actual determinant.
3. **`right_dn_orbs` crash (found during testing).** Its bounds guard reused
   the *up*-channel electron count (`n = up_orbs.size()`) instead of
   `dn_orbs.size()`. For a wavefunction with `n_up != n_dn` and `AS_last`
   falling between the two counts, this walked an iterator more than one
   past `dn_orbs.end()` (undefined behavior), observed as a
   `std::length_error` crash in testing.

Since `Compute_overlap` only ever consumes `Gs_index_map` (never
`AS_first`/`AS_last` directly), all three bugs corrupted both `Lambda` and
the projected overlap, even though `Compute_overlap` never uses `Lambda`'s
value for normalization.

### Fix

Replaced the position-based core/virtual split with a value-based one.
`AS_VBM_up` / `AS_VBM_dn` — the highest `A_S` orbital filled when `A_S`'s
own sorted orbital list is filled with `n_up_s` / `n_dn_s` electrons from
the bottom, i.e. `A_S`'s own HOMO in the reference state, per spin — now
define the boundary: any non-`A_S` orbital below its spin's VBM must be
occupied (reference core, gap orbitals included), any non-`A_S` orbital
above it must be empty.

This correctly handles non-contiguous `A_S` and removes all position-based
iterator arithmetic on `up_orbs`/`dn_orbs`, eliminating all three bugs
above. `e_trans` / `ne_s_up` / `ne_s_dn` (via `set_intersection`) were
already correct and are unchanged. `Compute_overlap`, `map_indices`, and
`wfc_proj.cc` were not touched — they only ever consume `Gs_index_map`, so
they inherit the fix automatically.

**Assumption relied on:** `A_S`/`A_L` band lists in the config are given in
ascending band order matching SHCI's own wavefunction orbital-index
convention (so `index_map` comes out already sorted).

**Known limitation:** `n_up_s == 0` (or `n_dn_s == 0`) — `A_S` empty for
that spin in the reference state — is not yet handled; `Compute_Lambda`
throws rather than silently producing a wrong answer.

**Verification:** compiled clean with `g++ -fsyntax-only -std=c++17`
against the real project headers. Empirically tested with a standalone
harness constructing `Wavefunction`/`Det` objects directly (no file I/O)
and calling `Compute_Lambda` on synthetic determinants: a clean
contiguous-`A_S` reference (correctly admitted), a leak-beyond-`AS_last`
case with a small electron count (correctly rejected — this is the bug that
crashed the original code under `n_up != n_dn`), a non-contiguous
`A_S = {2,5}` with gap orbitals correctly filled per the VBM rule
(correctly admitted), the same setup with a hole in a gap orbital — `A_S`'s
own count still exactly right (`e_trans == 0`) but a required-core gap
orbital left empty (correctly rejected — the "`e_trans == 0` but there's a
transition within `A_L − A_S`" case), and a gap-filled-but-with-an-extra-leak
variant (correctly rejected). All cases passed against the patched code; the
`n_up != n_dn` leak case crashed the original (pre-fix) code as described
above.
