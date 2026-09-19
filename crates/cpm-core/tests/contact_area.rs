use cpm_core::lattice::{Boundary, Lattice, Neighborhood};
use cpm_core::world::World;
use cpm_core::MEDIUM;
use std::collections::HashMap;

// 3x3 grid of single-pixel cells (each pixel its own cell id), Moore
// neighborhood (order=2 -> 8-connected in 2D, matching this project's
// influenza-sego2022 params.yaml `neighbor_order: 3`, which clamps to the
// same 8-connected set in 2D). The center cell is fully surrounded -> all 8
// Moore neighbors are OTHER cells' pixels, all initially type A(1).
//
// NOTE: `cell.surface` (recompute_trackers) counts unlike-owner pixels over
// `lattice.neighbors` (the full configured CPM neighborhood), not the
// axis-only `lattice.face_neighbors`. For the required invariant
// (sum of per-type contact areas == cell_surfaces()[cell_id]) to hold under
// this investigation's neighbor_order=3 config, `cell_contact_area_by_type`
// must use the same `lattice.neighbors` set that the surface tracker uses.
fn setup() -> (World, u32, [[u32; 3]; 3]) {
    let lat = Lattice::new([5, 5, 1], [Boundary::NoFlux; 3], Neighborhood::new(false, 2));
    let mut w = World::new(lat, 10.0);
    let mut ids = [[0u32; 3]; 3];
    for y in 0..3 {
        for x in 0..3 {
            let id = w.add_cell(1, 1.0, 1.0, 0.0, 0.0);
            let idx = w.lattice.index(x + 1, y + 1, 0);
            w.paint(idx, id);
            ids[y][x] = id;
        }
    }
    w.recompute_trackers();
    let center = ids[1][1];
    (w, center, ids)
}

#[test]
fn all_type_a_neighbors_tally_under_type_a_and_matches_surface() {
    let (w, center, _ids) = setup();
    let map = w.cell_contact_area_by_type(center);
    let mut expected: HashMap<u16, i64> = HashMap::new();
    expected.insert(1, 8);
    assert_eq!(map, expected);
    let sum: i64 = map.values().sum();
    assert_eq!(sum, w.cells[center as usize].surface);
    assert_eq!(sum, 8);
}

#[test]
fn flipping_one_neighbor_to_type_b_splits_the_tally_but_keeps_the_sum() {
    let (mut w, center, ids) = setup();
    // Corner neighbor at grid-local (0,0) -> retype to B(2).
    let corner = ids[0][0];
    w.set_cell_type(corner, 2);

    let map = w.cell_contact_area_by_type(center);
    let mut expected: HashMap<u16, i64> = HashMap::new();
    expected.insert(1, 7);
    expected.insert(2, 1);
    assert_eq!(map, expected);

    let sum: i64 = map.values().sum();
    assert_eq!(sum, w.cells[center as usize].surface);
    assert_eq!(sum, 8);
}

#[test]
fn medium_contacts_are_counted_under_type_medium() {
    // A single cell with no neighbors at all: every unlike-owner face is
    // medium (type 0). Use a lone 3x3-pixel cell in the middle of a NoFlux
    // 5x5 lattice with Moore neighborhood: hand-verified surface = 32 (see
    // world.rs::surface_of_isolated_3x3_moore_is_correct), all of it medium.
    let lat = Lattice::new([5, 5, 1], [Boundary::NoFlux; 3], Neighborhood::new(false, 2));
    let mut w = World::new(lat, 10.0);
    let a = w.add_cell(1, 9.0, 1.0, 0.0, 0.0);
    for y in 1..4 {
        for x in 1..4 {
            let idx = w.lattice.index(x, y, 0);
            w.paint(idx, a);
        }
    }
    w.recompute_trackers();

    let map = w.cell_contact_area_by_type(a);
    let mut expected: HashMap<u16, i64> = HashMap::new();
    expected.insert(MEDIUM as u16, 32);
    assert_eq!(map, expected);
    assert_eq!(map.values().sum::<i64>(), w.cells[a as usize].surface);
}
