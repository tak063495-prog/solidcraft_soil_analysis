use nalgebra::{DVector, Vector3};
use soil_fem_3d::{Result, joints, mesh::Mesh, rigidity, stages::StageCase};
use std::collections::BTreeSet;

fn fixture(hinge: bool) -> Result<(Mesh, joints::Configuration)> {
    let mesh = serde_json::from_str(if hinge {
        include_str!("../examples/joint-crack-tip-hinge-mesh.json")
    } else {
        include_str!("../examples/joint-crack-tip-mesh.json")
    })?;
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-crack-tip.json"))?;
    let mut cfg = case.joint_system.unwrap();
    if hinge {
        cfg.joints[0].crack_tip = Some(serde_json::from_str(include_str!(
            "../examples/joint-crack-tip-hinge-declaration.json"
        ))?);
    }
    Ok((mesh, cfg))
}

#[test]
fn declared_tip_cancels_exactly_and_preserves_virtual_work_force_and_moment() -> Result<()> {
    let (mesh, cfg) = fixture(false)?;
    let tip: BTreeSet<_> = cfg.joints[0]
        .crack_tip
        .as_ref()
        .unwrap()
        .edges
        .iter()
        .flatten()
        .copied()
        .collect();
    let model = joints::Model::new(&mesh, cfg)?;
    let layout = model.layout();
    let u: Vec<[f64; 3]> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, p)| {
            [
                1e-6 * (i as f64).sin(),
                2e-6 * (i as f64).cos(),
                -1e-6 * p[0] + 1e-6 * (i as f64).sin(),
            ]
        })
        .collect();
    let virtual_u: Vec<[f64; 3]> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, p)| {
            [
                0.1 + (i as f64).sin(),
                -0.2 + p[0] * p[1],
                p[2] + (i as f64).cos(),
            ]
        })
        .collect();
    let active = vec![true; mesh.tetrahedra.len()];
    let eval = model.evaluate(&model.initial_state()?, &u, &active, None, true)?;
    let mut expected_work = 0.;
    for ((p, s), coefficients) in layout
        .iter()
        .zip(&eval.state.points)
        .zip(model.jump_coefficients())
    {
        let mut jump = Vector3::zeros();
        let mut virtual_jump = Vector3::zeros();
        // Independently evaluate the Tri6 polynomial at the GP natural point,
        // recovered from its affine world position (straight fixture only).
        let a = Vector3::from(mesh.nodes[p.nodes_a[0]]);
        let e = nalgebra::Matrix3x2::from_columns(&[
            Vector3::from(mesh.nodes[p.nodes_a[1]]) - a,
            Vector3::from(mesh.nodes[p.nodes_a[2]]) - a,
        ]);
        let r = (e.transpose() * e).try_inverse().unwrap()
            * e.transpose()
            * (Vector3::from(p.position_m) - a);
        let l = [1. - r[0] - r[1], r[0], r[1]];
        let shape = [
            l[0] * (2. * l[0] - 1.),
            l[1] * (2. * l[1] - 1.),
            l[2] * (2. * l[2] - 1.),
            4. * l[0] * l[1],
            4. * l[1] * l[2],
            4. * l[2] * l[0],
        ];
        for (i, &w) in shape.iter().enumerate() {
            jump += w * (Vector3::from(u[p.nodes_b[i]]) - Vector3::from(u[p.nodes_a[i]]));
            virtual_jump += w
                * (Vector3::from(virtual_u[p.nodes_b[i]]) - Vector3::from(virtual_u[p.nodes_a[i]]));
            assert!((w - p.shape[i]).abs() < 1e-12);
        }
        assert!((jump - Vector3::from(s.jump_m)).norm() < 1e-15);
        expected_work += p.weight_m2 * virtual_jump.dot(&Vector3::from(s.traction_kpa));
        assert!(coefficients.iter().all(|(n, _)| !tip.contains(n)));
    }
    let work: f64 = eval
        .internal_kn
        .iter()
        .zip(virtual_u.iter().flatten())
        .map(|(f, u)| f * u)
        .sum();
    assert!((work - expected_work).abs() < 1e-12);
    assert!(
        model
            .connectivity()
            .iter()
            .flatten()
            .all(|g| !tip.contains(&(g / 3)))
    );
    for n in tip {
        assert!(eval.internal_kn[3 * n..3 * n + 3].iter().all(|f| *f == 0.));
        assert!(
            eval.tangent
                .iter()
                .all(|(g, h, _)| g / 3 != n && h / 3 != n)
        );
    }
    assert!(Vector3::from(eval.report.resultant_internal_force_kn).norm() < 1e-10);
    assert!(Vector3::from(eval.report.resultant_internal_moment_kn_m).norm() < 1e-10);
    model.validate_state(&eval.state, &u, &active)?;
    Ok(())
}

#[test]
fn compression_opening_and_reversed_shear_have_consistent_tangents() -> Result<()> {
    let (mesh, cfg) = fixture(true)?;
    let upper: BTreeSet<_> = mesh
        .region_ids
        .iter()
        .enumerate()
        .filter(|(_, r)| **r == 2)
        .flat_map(|(e, _)| mesh.cell_nodes(e))
        .collect();
    let lower: BTreeSet<_> = mesh
        .region_ids
        .iter()
        .enumerate()
        .filter(|(_, r)| **r == 1)
        .flat_map(|(e, _)| mesh.cell_nodes(e))
        .collect();
    let model = joints::Model::new(&mesh, cfg)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let initial = model.initial_state()?;
    for (shear, opening) in [(1e-5, -1e-4), (0., 1e-4), (-1e-3, -1e-4)] {
        let u: Vec<_> = mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(n, p)| {
                if upper.contains(&n) && !lower.contains(&n) {
                    [shear * p[0], 0., opening * p[0]]
                } else {
                    [0.; 3]
                }
            })
            .collect();
        let eval = model.evaluate(&initial, &u, &active, None, true)?;
        let direction: Vec<[f64; 3]> = (0..u.len())
            .map(|n| [(n as f64).sin(), (n as f64).cos(), 0.1])
            .collect();
        let h = 1e-9;
        let shifted = |sign: f64| {
            u.iter()
                .zip(&direction)
                .map(|(v, d)| std::array::from_fn(|a| v[a] + sign * h * d[a]))
                .collect::<Vec<_>>()
        };
        let plus = model.evaluate(&initial, &shifted(1.), &active, None, false)?;
        let minus = model.evaluate(&initial, &shifted(-1.), &active, None, false)?;
        let fd = DVector::from_iterator(
            3 * u.len(),
            plus.internal_kn
                .iter()
                .zip(&minus.internal_kn)
                .map(|(a, b)| (a - b) / (2. * h)),
        );
        let mut predicted = DVector::zeros(fd.len());
        for (g, k, v) in eval.tangent {
            predicted[g] += v * direction[k / 3][k % 3];
        }
        if fd.norm() > 0. {
            assert!((&predicted - &fd).norm() / fd.norm() < 1e-7);
        } else {
            assert!(predicted.norm() < 1e-10);
        }
    }
    Ok(())
}

#[test]
fn shared_tip_edge_keeps_one_open_hinge_and_contact_only_closes_it_in_compression() -> Result<()> {
    let (mesh, cfg) = fixture(true)?;
    let lower: BTreeSet<_> = mesh
        .region_ids
        .iter()
        .enumerate()
        .filter(|(_, r)| **r == 1)
        .flat_map(|(e, _)| mesh.cell_nodes(e))
        .collect();
    let active = vec![true; mesh.tetrahedra.len()];
    let geometry = rigidity::Geometry::new(&mesh, &active)?;
    let mut next = 0;
    let rows: Vec<_> = (0..3 * mesh.nodes.len())
        .map(|g| {
            if lower.contains(&(g / 3)) {
                vec![]
            } else {
                let v = vec![(next, 1.)];
                next += 1;
                v
            }
        })
        .collect();
    let model = joints::Model::new(&mesh, cfg)?;
    for (theta, missing) in [(0., 1), (1e-4, 1), (-1e-4, 0)] {
        let u: Vec<_> = mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(n, p)| {
                if lower.contains(&n) {
                    [0.; 3]
                } else {
                    [-theta * (p[2] - 0.5), 0., theta * p[0]]
                }
            })
            .collect();
        let eval = model.evaluate(&model.initial_state()?, &u, &active, None, false)?;
        let mut space = geometry.constrained(&rows)?;
        assert_eq!(space.report().continuum_components, 2);
        assert_eq!(space.report().unrestrained_modes, 1);
        model.add_rigid_support(&eval.state, &active, 1e-8, &mut space)?;
        assert_eq!(space.report().unrestrained_modes, missing);
        assert_eq!(space.require_supported().is_ok(), missing == 0);
    }
    Ok(())
}

#[test]
fn undeclared_partial_wrong_owner_midpoint_duplicate_and_fully_tied_faces_are_rejected()
-> Result<()> {
    let (mesh, cfg) = fixture(false)?;
    let mut bad = cfg.clone();
    bad.joints[0].crack_tip = None;
    assert!(joints::Model::new(&mesh, bad).is_err());
    for mode in 0..6 {
        let mut bad = cfg.clone();
        let tip = bad.joints[0].crack_tip.as_mut().unwrap();
        match mode {
            0 => tip.edges[0][2] = tip.edges[0][0],
            1 => tip.edges.push(tip.edges[0]),
            2 => tip.faces[0].owner_a = tip.faces[0].owner_b,
            3 => {
                tip.faces.pop();
            }
            4 => tip.faces[0].face_b = tip.faces[0].face_a,
            _ => tip.edges[0].swap(1, 2),
        }
        assert!(joints::Model::new(&mesh, bad).is_err(), "mode {mode}");
    }
    let mut bad = cfg;
    bad.joints[0].side_b = bad.joints[0].side_a.clone();
    assert!(joints::Model::new(&mesh, bad).is_err());
    Ok(())
}

#[test]
fn rehashed_checkpoint_cannot_change_tip_surface_or_owning_sides() -> Result<()> {
    use soil_fem_3d::{
        checkpoint::{Checkpoint, hash_json},
        ground::Ground,
        stages,
    };
    let (mesh, _) = fixture(false)?;
    let ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-joint-crack-tip.json"))?;
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-crack-tip.json"))?;
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let saved = serde_json::to_vec(&Checkpoint::seal(&mesh, partial)?)?;
    for mode in 0..4 {
        let mut bad: Checkpoint = serde_json::from_slice(&saved)?;
        let tip = bad.run.case.joint_system.as_mut().unwrap().joints[0]
            .crack_tip
            .as_mut()
            .unwrap();
        match mode {
            0 => {
                tip.edges.pop();
            }
            1 => {
                tip.faces.pop();
            }
            2 => {
                let f = &mut tip.faces[0];
                std::mem::swap(&mut f.owner_a, &mut f.owner_b);
            }
            _ => {
                bad.run.joint_integration_points[0].owner_a += 1;
            }
        }
        bad.run_sha256 = hash_json(&bad.run)?;
        assert!(
            bad.resume(&mesh, ground.clone(), case.clone(), None)
                .is_err(),
            "mode={mode}"
        );
    }
    Ok(())
}

#[test]
fn declared_tip_rejects_geometric_ambiguity_even_after_an_opposing_face_was_used() -> Result<()> {
    let tolerance = 1e-10;
    let mut nodes = Vec::new();
    let mut tetrahedra = Vec::new();
    let mut faces_a = Vec::new();
    let mut faces_b = Vec::new();
    for shift in [0., 0.9 * tolerance] {
        let n = nodes.len();
        nodes.extend([
            [shift, 0., 0.],
            [1. + shift, 0., 0.],
            [shift, 1., 0.],
            [shift, 0., -1.],
            [shift + 0.6 * tolerance, 1., 0.],
            [shift, 0., 1.],
        ]);
        tetrahedra.extend([[n, n + 2, n + 1, n + 3], [n, n + 1, n + 4, n + 5]]);
        faces_a.push([n, n + 1, n + 2]);
        faces_b.push([n, n + 1, n + 4]);
    }
    let linear: Mesh = serde_json::from_value(serde_json::json!({
        "schema":"soil-fem-mesh/2","units":"m-kN-kPa","axis":"z_up","crs":"local","origin_m":[0.,0.,0.],
        "nodes":nodes,"tetrahedra":tetrahedra,"material_ids":[1,1,1,1],"region_ids":[1,1,1,1],
        "node_sets":{},"face_sets":{"Seam:a":faces_a,"Seam:b":faces_b},
        "regions":[{"id":1,"body":"same","material_id":1,"bounds_m":[[0.,0.,-1.],[1.1,1.,1.]],"volume_m3":2./3.,"source_sha256":"a".repeat(64)}]
    }))?;
    let mesh = linear.elevate_straight()?;
    let (_, mut cfg) = fixture(false)?;
    cfg.geometric_tolerance_m = tolerance;
    let mut edges = Vec::new();
    let mut faces = Vec::new();
    for k in 0..2 {
        let n = 6 * k;
        let cell = mesh.tetrahedra10[2 * k];
        let edge = soil_fem_3d::elements::EDGES
            .iter()
            .position(|&[i, j]| BTreeSet::from([cell[i], cell[j]]) == BTreeSet::from([n, n + 1]))
            .unwrap();
        edges.push([n, n + 1, cell[edge + 4]]);
        faces.push(joints::OwnedFace {
            face_a: faces_a[k],
            face_b: faces_b[k],
            owner_a: 2 * k,
            owner_b: 2 * k + 1,
        });
    }
    cfg.joints[0].crack_tip = Some(joints::CrackTip { edges, faces });
    // A1 matches only B1. A2 matches both B1 and B2. Reusing the first
    // assignment to hide that ambiguity must never make the input valid.
    assert!(joints::Model::new(&mesh, cfg).is_err());
    Ok(())
}

#[test]
fn a_branched_front_cannot_be_hidden_in_three_separate_joint_definitions() -> Result<()> {
    let mut nodes = vec![[0., 0., 0.]];
    let mut tets = Vec::new();
    let mut sets = serde_json::Map::new();
    let arms = [[1., 0.], [0., 1.], [-1., 0.]];
    for (k, [x, y]) in arms.into_iter().enumerate() {
        let n = nodes.len();
        nodes.extend([
            [x, y, 0.],
            [-y, x, 0.],
            [0., 0., -1.],
            [-y, x, 0.],
            [0., 0., 1.],
        ]);
        tets.extend([[0, n + 1, n, n + 2], [0, n, n + 3, n + 4]]);
        sets.insert(format!("J{k}:a"), serde_json::json!([[0, n, n + 1]]));
        sets.insert(format!("J{k}:b"), serde_json::json!([[0, n, n + 3]]));
    }
    let linear: Mesh = serde_json::from_value(serde_json::json!({
        "schema":"soil-fem-mesh/2","units":"m-kN-kPa","axis":"z_up","crs":"local","origin_m":[0.,0.,0.],
        "nodes":nodes,"tetrahedra":tets,"material_ids":[1,1,1,1,1,1],"region_ids":[1,1,1,1,1,1],
        "node_sets":{},"face_sets":sets,"regions":[{"id":1,"body":"same","material_id":1,"bounds_m":[[-1.,-1.,-1.],[1.,1.,1.]],"volume_m3":1.,"source_sha256":"a".repeat(64)}]
    }))?;
    let mesh = linear.elevate_straight()?;
    let (_, mut cfg) = fixture(false)?;
    let template = cfg.joints[0].clone();
    cfg.joints.clear();
    for k in 0..3 {
        let n = 1 + 5 * k;
        let cell = mesh.tetrahedra10[2 * k];
        let ei = soil_fem_3d::elements::EDGES
            .iter()
            .position(|&[i, j]| BTreeSet::from([cell[i], cell[j]]) == BTreeSet::from([0, n]))
            .unwrap();
        let mut d = template.clone();
        d.id = format!("J{k}");
        d.side_a = format!("J{k}:a");
        d.side_b = format!("J{k}:b");
        d.crack_tip = Some(joints::CrackTip {
            edges: vec![[0, n, cell[ei + 4]]],
            faces: vec![joints::OwnedFace {
                face_a: [0, n, n + 1],
                face_b: [0, n, n + 3],
                owner_a: 2 * k,
                owner_b: 2 * k + 1,
            }],
        });
        cfg.joints.push(d);
    }
    assert!(joints::Model::new(&mesh, cfg).is_err());
    Ok(())
}
