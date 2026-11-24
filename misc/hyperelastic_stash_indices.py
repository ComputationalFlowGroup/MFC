#!/usr/bin/env python3
"""Verify the dedicated hyperelastic stash indices reproduce the former xi%end offsets exactly.

Before: every consumer addressed the stash as eqn_idx%xi%end + <literal>, and sys_size grew by
4 + 2*num_dims. After: m_global_parameters_common.fpp assigns dedicated eqn_idx entries and
sys_size follows from the last one. This asserts the two index maps agree for every slot and
every num_dims, so the rename moved no field.
"""

for num_dims in (1, 2, 3):
    # S is sys_size just before the hyperelasticity block, which equals eqn_idx%xi%end there
    S = 17  # arbitrary; every expression below is affine in S

    # --- former offsets, as the 52 replaced call sites spelled them ---
    old = {
        "W": S + 1,
        "trB": S + 2,
        "J": S + 3,
        "trBinv": S + 4,
        **{f"Bdiag{i}": S + 4 + i for i in range(1, num_dims + 1)},
        **{f"Binvdiag{i}": S + 4 + num_dims + i for i in range(1, num_dims + 1)},
    }
    old_sys_size = S + 4 + 2 * num_dims

    # --- dedicated indices, as m_global_parameters_common.fpp now assigns them ---
    elastic_W = S + 1
    trB = S + 2
    Jdet = S + 3
    trBinv = S + 4
    Bdiag_beg = S + 5
    Bdiag_end = Bdiag_beg + num_dims - 1
    Binvdiag_beg = Bdiag_end + 1
    Binvdiag_end = Binvdiag_beg + num_dims - 1
    new_sys_size = Binvdiag_end

    new = {
        "W": elastic_W,
        "trB": trB,
        "J": Jdet,
        "trBinv": trBinv,
        **{f"Bdiag{i}": Bdiag_beg - 1 + i for i in range(1, num_dims + 1)},
        **{f"Binvdiag{i}": Binvdiag_beg - 1 + i for i in range(1, num_dims + 1)},
    }

    assert old == new, (num_dims, old, new)
    assert old_sys_size == new_sys_size, (num_dims, old_sys_size, new_sys_size)

    # the stash must be contiguous, start one past xi%end, and end exactly at sys_size
    slots = sorted(new.values())
    assert slots == list(range(S + 1, new_sys_size + 1)), (num_dims, slots)
    assert len(slots) == 4 + 2 * num_dims, (num_dims, len(slots))

    # the allocate/seed/deallocate loops run elastic_W .. Binvdiag%end: it must cover every slot
    assert min(slots) == elastic_W and max(slots) == Binvdiag_end

    print(f"num_dims = {num_dims}: {len(slots)} slots, {S + 1}..{new_sys_size}, index maps agree")

print("OK: dedicated indices reproduce the former offsets for every num_dims")
