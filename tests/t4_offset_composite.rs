use nalgebra::{DVector, Vector3};
use soil_fem_3d::{
    Result, elements, joint_law,
    mesh::{Mesh, Region},
    offset_composite::{
        Contact, Face, Patch, Plate,
        path::{self, Case, Placement, Soil, Step},
    },
    plates::Section,
};
use std::collections::BTreeMap;

fn fixture() -> Result<(Mesh, Case)> {
    // Upper face 0,1,2 of lower tet and matching lower face 4,5,6 of upper tet.
    let regions = (1..=2)
        .map(|id| Region {
            id,
            body: format!("soil{id}"),
            material_id: 1,
            bounds_m: [
                [0., 0., if id == 1 { -1. } else { 0.1 }],
                [1., 1., if id == 1 { 0. } else { 1.1 }],
            ],
            volume_m3: 1. / 6.,
            source_sha256: "0".repeat(64),
        })
        .collect();
    let mesh = Mesh {
        schema: "soil-fem-mesh/2".into(),
        units: "m-kN-kPa".into(),
        axis: "z_up".into(),
        crs: "local".into(),
        origin_m: [0.; 3],
        nodes: vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., -1.],
            [0., 0., 0.1],
            [1., 0., 0.1],
            [0., 1., 0.1],
            [0., 0., 1.1],
        ],
        tetrahedra: vec![[0, 2, 1, 3], [4, 5, 6, 7]],
        tetrahedra10: vec![],
        quadrature_order: None,
        material_ids: vec![1, 1],
        node_sets: BTreeMap::new(),
        face_sets: BTreeMap::new(),
        regions,
        region_ids: vec![1, 2],
        excavation_sets: BTreeMap::new(),
    }
    .elevate_straight()?;
    let face = |cell: usize, ids: [usize; 3]| {
        let tet = mesh.tetrahedra10[cell];
        let mut nodes = [0; 6];
        nodes[..3].copy_from_slice(&ids);
        for (k, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
            let edge = elements::EDGES
                .iter()
                .position(|[i, j]| {
                    [tet[*i], tet[*j]].contains(&ids[a]) && [tet[*i], tet[*j]].contains(&ids[b])
                })
                .unwrap();
            nodes[k + 3] = tet[4 + edge];
        }
        Face { cell, nodes }
    };
    let contact = Contact {
        id: "exposed-B".into(),
        soil_b: face(1, [4, 5, 6]),
        material: joint_law::Material {
            normal_stiffness_kpa_m: 1000.,
            tangential_stiffness_kpa_m: 500.,
            friction_coefficient: 0.3,
        },
    };
    let plate = Plate {
        id: "liner".into(),
        section: Section {
            id: "P".into(),
            young_kpa: 200000.,
            poisson: 0.25,
            thickness_m: 0.1,
            unit_weight_kn_m3: 0.,
            drilling_factor: 0.001,
        },
        soil_a: face(0, [0, 1, 2]),
        offset_m: -0.05,
        exposed_side: 1,
        tie_at_physical_surface: true,
        contacts: vec![contact],
    };
    Ok((
        mesh,
        Case {
            plate,
            soil: BTreeMap::from([(
                1,
                Soil {
                    young_kpa: 10000.,
                    poisson: 0.25,
                },
            )]),
            steps: vec![],
        },
    ))
}
fn place(p: Placement) -> Step {
    Step::Install {
        contacts: BTreeMap::from([("exposed-B".into(), p)]),
        holder: "displacement supports".into(),
    }
}
fn load(d: Vec<f64>) -> Step {
    Step::Load {
        free_dofs: vec![],
        nodal_load: vec![0.; d.len()],
        prescribed_configuration: d,
        plate_traction_kpa: [0.; 3],
    }
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t, "{a} vs {b}");
}

#[test]
fn face_ownership_and_open_close() -> Result<()> {
    let (mesh, mut case) = fixture()?;
    let plate = &case.plate;
    let mut duplicate = plate.clone();
    duplicate.contacts.push(plate.contacts[0].clone());
    duplicate.contacts[1].id = "redeclared".into();
    assert!(Patch::new(&mesh, &duplicate).is_err());
    duplicate.contacts.reverse();
    assert!(Patch::new(&mesh, &duplicate).is_err());
    let mut same = plate.clone();
    same.contacts[0].soil_b = same.soil_a.clone();
    assert!(Patch::new(&mesh, &same).is_err());
    let mut missing = plate.clone();
    missing.contacts[0].soil_b.cell = 8;
    assert!(Patch::new(&mesh, &missing).is_err());
    let mut steps = vec![place(Placement::StressFree)];
    for (dz, dx) in [
        (-0.01, 0.001),
        (0.002, 0.001),
        (-0.01, 0.001),
        (-0.01, 0.02),
    ] {
        let mut d = vec![0.; mesh.nodes.len() * 3 + 18];
        for n in mesh.tetrahedra10[1] {
            d[3 * n + 2] = dz;
            d[3 * n] = dx;
        }
        steps.push(load(d));
    }
    case.steps = steps;
    let run = path::run(&mesh, &case, 5, None)?;
    for p in &run.frames[5].active.as_ref().unwrap().contacts[0].points {
        assert_eq!(p.status, joint_law::Status::Slip);
        close(p.traction_kpa[0], 3., 1e-12);
        assert!(p.plastic_dissipation_kpa_m > 0.);
    }
    for (i, q) in [(2, -10.), (3, 0.), (4, -10.)] {
        let f = &run.frames[i];
        for p in &f.active.as_ref().unwrap().contacts[0].points {
            close(p.traction_kpa[2], q, 1e-11);
        }
        for n in case.plate.soil_a.nodes {
            assert_eq!(f.configuration[3 * n], 0.);
        }
        let resultant: f64 = case.plate.contacts[0]
            .soil_b
            .nodes
            .iter()
            .map(|n| f.reactions[3 * n + 2])
            .sum();
        close(resultant, 0.5 * q, 1e-10);
    }
    Ok(())
}

#[test]
fn signed_arms_reversal_and_rigid_rotation() -> Result<()> {
    let (mesh, case) = fixture()?;
    let mut plate = case.plate;
    let patch = Patch::new(&mesh, &plate)?;
    let theta = Vector3::new(0.01, -0.02, 0.003);
    let shift = Vector3::new(0.002, 0.005, -0.001);
    let mut d = DVector::zeros(54);
    for i in 0..6 {
        let x = Vector3::from(mesh.nodes[plate.soil_a.nodes[i]]);
        d.rows_mut(6 * i, 3).copy_from(&(shift + theta.cross(&x)));
        d.rows_mut(6 * i + 3, 3).copy_from(&theta);
        d.rows_mut(36 + 3 * i, 3)
            .copy_from(&(shift + theta.cross(&(x + patch.arm))));
    }
    for p in &patch.points {
        assert!((&p.jump * &d).norm() < 1e-16);
    }
    let section = patch.plate.recover("p".into(), &d.rows(0, 36).into_owned());
    for s in section.samples {
        assert!(
            s.membrane_strain
                .iter()
                .chain(&s.curvature_per_m)
                .chain(&s.transverse_shear_strain)
                .all(|v| v.abs() < 1e-15)
        );
    }
    for f in [&mut plate.soil_a, &mut plate.contacts[0].soil_b] {
        f.nodes.swap(1, 2);
        f.nodes.swap(3, 5);
    }
    plate.offset_m = 0.05;
    plate.exposed_side = -1;
    let reversed = Patch::new(&mesh, &plate)?;
    assert_eq!(patch.arm, reversed.arm);
    plate.exposed_side = 1;
    assert!(Patch::new(&mesh, &plate).is_err());
    Ok(())
}

#[test]
fn transfer_checks_every_arm_moment_and_freezes_generations() -> Result<()> {
    let (mesh, mut case) = fixture()?;
    let patch = Patch::new(&mesh, &case.plate)?;
    // For constant traction, corner integrals vanish. Each midside has area/3.
    let q = Vector3::new(1., 0., -10.);
    let mut actions = vec![0.; 54];
    for i in 3..6 {
        for a in 0..3 {
            actions[6 * i + a] = -q[a] / 6.;
            actions[6 * i + 3 + a] = -patch.arm.cross(&q)[a] / 6.;
            actions[36 + 3 * i + a] = q[a] / 6.;
        }
    }
    let transfer = |q: Vector3<f64>, actions: Vec<f64>| {
        place(Placement::TransferForce {
            supplier: "jack".into(),
            traction_kpa: vec![q.into(); patch.points.len()],
            requested_nodal_actions: actions,
        })
    };
    let mut d = vec![0.; 3 * mesh.nodes.len() + 18];
    d[0] = 0.003;
    d[mesh.nodes.len() * 3 + 1] = 0.01;
    let mut bad = actions.clone();
    bad[21 + 1] += 0.001;
    case.steps = vec![load(d), transfer(q, bad)];
    assert!(path::run(&mesh, &case, 2, None).is_err());
    for inadmissible in [Vector3::new(0., 0., 1.), Vector3::new(4., 0., -10.)] {
        case.steps[1] = transfer(inadmissible, actions.clone());
        assert!(path::run(&mesh, &case, 2, None).is_err());
    }
    let mut cancellation = actions.clone();
    cancellation[0] += 0.1;
    cancellation[6] -= 0.1;
    case.steps[1] = transfer(q, cancellation);
    assert!(path::run(&mesh, &case, 2, None).is_err());
    case.steps[1] = transfer(q, actions);
    case.steps.push(Step::Remove {
        holder: "shoring".into(),
    });
    case.steps.push(place(Placement::StressFree));
    let run = path::run(&mesh, &case, 4, None)?;
    let g = run.frames[2].active.as_ref().unwrap();
    assert_eq!(g.contacts[0].endpoint_work_kn_m, 0.);
    close(
        g.last_sections.as_ref().unwrap().bending_energy_kn_m,
        0.,
        1e-20,
    );
    assert_eq!(run.frames[4].active.as_ref().unwrap().generation, 2);
    assert_eq!(
        serde_json::to_vec(g)?,
        serde_json::to_vec(&run.frames[4].archived[0])?
    );
    for i in 2..=4 {
        let mut cp = path::run(&mesh, &case, i, None)?;
        let resumed = path::run(&mesh, &case, 4, Some(&cp))?;
        assert_eq!(serde_json::to_vec(&run)?, serde_json::to_vec(&resumed)?);
        if i >= 3 {
            cp.frames[i].archived[0].installation[0] += 0.01;
            assert!(path::run(&mesh, &case, 4, Some(&cp)).is_err());
        }
    }
    Ok(())
}

#[test]
fn instantaneous_and_accepted_work_use_different_origins() -> Result<()> {
    let (mesh, mut case) = fixture()?;
    let mut d = vec![0.; mesh.nodes.len() * 3 + 18];
    case.steps.push(place(Placement::StressFree));
    for dz in [-0.01, 0., -0.01] {
        for n in mesh.tetrahedra10[1] {
            d[3 * n + 2] = dz;
        }
        case.steps.push(load(d.clone()));
    }
    let run = path::run(&mesh, &case, 4, None)?;
    for i in [2, 4] {
        let f = &run.frames[i];
        close(f.instantaneous_external_kn_m, 0.05, 1e-12);
        close(f.cumulative_reaction_work_kn_m, 0.025, 1e-12);
    }
    close(run.frames[3].cumulative_reaction_work_kn_m, 0., 1e-12);
    let mut replay = path::run(&mesh, &case, 2, None)?;
    replay.frames[2].cumulative_reaction_work_kn_m *= 2.;
    assert!(path::run(&mesh, &case, 4, Some(&replay)).is_err());
    Ok(())
}

#[test]
fn point_work_all_nodal_forces_and_tangent_are_independent() -> Result<()> {
    let (mesh, case) = fixture()?;
    let patch = Patch::new(&mesh, &case.plate)?;
    let material = &case.plate.contacts[0].material;
    let old: Vec<_> = patch
        .points
        .iter()
        .map(|_| material.birth(patch.normal.into(), [0.; 3], [0.; 3]))
        .collect::<Result<_>>()?;
    let mut d = DVector::zeros(54);
    for i in 0..6 {
        d[36 + 3 * i] = 0.001;
        d[36 + 3 * i + 2] = -0.01;
    }
    let trial = patch.assemble_contact(material, &old, &d)?;
    for col in 0..54 {
        let mut probe = DVector::zeros(54);
        probe[col] = 1.;
        let mut local_work = 0.;
        let mut force = DVector::zeros(54);
        let mut tangent = DVector::zeros(54);
        for (p, s) in patch.points.iter().zip(&trial.points) {
            let mut du = Vector3::zeros();
            for (i, n) in p.shape.iter().enumerate() {
                let v = |k| Vector3::new(probe[k], probe[k + 1], probe[k + 2]);
                du += (v(36 + 3 * i) - v(6 * i) - v(6 * i + 3).cross(&patch.arm)) * *n;
            }
            let q = Vector3::from(s.traction_kpa);
            let dq = material.tangential_stiffness_kpa_m * du
                + (material.normal_stiffness_kpa_m - material.tangential_stiffness_kpa_m)
                    * patch.normal
                    * du.dot(&patch.normal);
            local_work += q.dot(&du) * p.weight;
            for (v, load) in [(&mut force, q), (&mut tangent, dq)] {
                for (i, n) in p.shape.iter().enumerate() {
                    let f = load * (p.weight * n);
                    for a in 0..3 {
                        v[6 * i + a] -= f[a];
                        v[6 * i + 3 + a] -= patch.arm.cross(&f)[a];
                        v[36 + 3 * i + a] += f[a];
                    }
                }
            }
        }
        close(local_work, trial.force.dot(&probe), 1e-14);
        assert!((&trial.force - force).amax() < 1e-12);
        assert!((&trial.tangent * &probe - tangent).amax() < 1e-8);
    }
    Ok(())
}

#[test]
fn ten_kn_neutral_plane_load_transports_its_moment_and_work() -> Result<()> {
    let (mesh, case) = fixture()?;
    let mut plate = case.plate;
    plate.contacts.clear();
    plate.tie_at_physical_surface = false;
    for e in [-0.05, 0., 0.05] {
        plate.offset_m = e;
        let patch = Patch::new(&mesh, &plate)?;
        let load = patch.plate.uniform_load([20., 0., 0.]);
        let mut total = DVector::zeros(6);
        let mut d = DVector::zeros(36);
        for i in 0..6 {
            total += load.rows(6 * i, 6);
            d[6 * i + 4] = 0.01;
        }
        close(total[0], 10., 1e-12);
        close(total[4], -e * 10., 1e-12);
        close(load.dot(&d), -e * 0.1, 1e-14);
        let perpendicular = patch.plate.uniform_load([0., 0., 20.]);
        close(perpendicular.dot(&d), 0., 1e-14);
    }
    Ok(())
}
