use cpm_core::lattice::{Boundary, Lattice, Neighborhood};
use cpm_core::world::World;

// One secreting cell (3x3 block), a field with a per-type rate, no diffusion/decay
// so secretion is the only source of mass. Compare total field concentration after
// N steps under default scale (1.0, unset), scale=0.5, and scale=0.0.
fn total_conc_after_steps(scale: Option<f64>, steps: u32) -> f64 {
    let lat = Lattice::new([6, 6, 1], [Boundary::NoFlux; 3], Neighborhood::new(false, 2));
    let mut w = World::new(lat, 10.0);
    let a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0);
    for y in 1..4 {
        for x in 1..4 {
            let idx = w.lattice.index(x, y, 0);
            w.paint(idx, a);
        }
    }
    w.recompute_trackers();
    let fi = w.add_field("virus", 0.0, 0.0);
    w.set_secretion(fi, 1, 10.0);
    if let Some(s) = scale {
        w.set_cell_secretion_scale(fi, a, s);
    }

    for _ in 0..steps {
        w.advance_fields();
    }

    w.fields[fi].conc.iter().map(|&c| c as f64).sum()
}

#[test]
fn default_scale_is_unset_and_numerically_identical_to_no_scale() {
    // Backward compatibility: not calling set_cell_secretion_scale at all must give
    // exactly the same result as an explicit scale of 1.0.
    let unset = total_conc_after_steps(None, 5);
    let explicit_one = total_conc_after_steps(Some(1.0), 5);
    assert_eq!(unset, explicit_one, "unset scale must equal explicit scale=1.0");
    assert!(unset > 0.0, "sanity: some secretion should have happened, got {unset}");
}

#[test]
fn scale_half_halves_total_secretion() {
    let full = total_conc_after_steps(None, 5);
    let half = total_conc_after_steps(Some(0.5), 5);
    assert!(
        (half - 0.5 * full).abs() < 1e-6,
        "half-scale secretion should be ~0.5x full: full={full}, half={half}"
    );
}

#[test]
fn scale_zero_suppresses_secretion() {
    let zero = total_conc_after_steps(Some(0.0), 5);
    assert!(zero.abs() < 1e-9, "scale=0.0 should give ~0 secretion, got {zero}");
}

#[test]
fn unrelated_cell_is_unaffected_by_another_cells_scale() {
    // Two secreting cells of the same type in the same field; scaling one down must
    // not affect the other's contribution.
    let lat = Lattice::new([10, 6, 1], [Boundary::NoFlux; 3], Neighborhood::new(false, 2));
    let mut w = World::new(lat, 10.0);
    let a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0);
    let b = w.add_cell(1, 9.0, 1.0, 0.0, 0.0);
    for y in 1..4 {
        for x in 1..4 {
            w.paint(w.lattice.index(x, y, 0), a);
        }
        for x in 6..9 {
            w.paint(w.lattice.index(x, y, 0), b);
        }
    }
    w.recompute_trackers();
    let fi = w.add_field("virus", 0.0, 0.0);
    w.set_secretion(fi, 1, 10.0);
    w.set_cell_secretion_scale(fi, a, 0.0);

    w.advance_fields();

    let a_pixel = w.lattice.index(2, 2, 0);
    let b_pixel = w.lattice.index(7, 2, 0);
    assert!(
        (w.fields[fi].conc[a_pixel] as f64).abs() < 1e-9,
        "a's secretion should be suppressed, got {}",
        w.fields[fi].conc[a_pixel]
    );
    assert!(
        (w.fields[fi].conc[b_pixel] as f64 - 10.0).abs() < 1e-6,
        "b's secretion should be untouched at 10.0, got {}",
        w.fields[fi].conc[b_pixel]
    );
}
