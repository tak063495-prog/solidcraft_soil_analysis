use nalgebra::{Matrix3, Vector3};
use serde_json::json;
use soil_fem_3d::{
    Result, beams,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    plates::{Kernel, Section, V18},
    stages::{self, StageCase},
};

fn section(thickness: f64, drilling: f64) -> Section {
    Section {
        id: "wall".into(),
        young_kpa: 2e8,
        poisson: 0.,
        thickness_m: thickness,
        unit_weight_kn_m3: 0.,
        drilling_factor: drilling,
    }
}

#[test]
fn dkt_membrane_curvature_patch_six_rigid_modes_and_energy_are_objective() -> Result<()> {
    let a = 0.43_f64;
    let rotation = Matrix3::new(a.cos(), -a.sin(), 0., a.sin(), a.cos(), 0., 0., 0., 1.)
        * Matrix3::new(1., 0., 0., 0., a.cos(), -a.sin(), 0., a.sin(), a.cos());
    let points = [[0., 0., 0.], [2., 0., 0.], [0.4, 1.3, 0.]];
    let kernel = Kernel::new(points, &section(0.02, 1e-4))?;
    let rotated_points =
        points.map(|x| (rotation * Vector3::from(x) + Vector3::new(2., -3., 4.)).into());
    let rotated = Kernel::new(rotated_points, &section(0.02, 1e-4))?;
    for mode in 0..6 {
        let rigid = V18::from_fn(|i, _| {
            let n = i / 6;
            let part = i % 6;
            if mode < 3 {
                if part == mode { 1. } else { 0. }
            } else {
                let mut w = Vector3::zeros();
                w[mode - 3] = 1.;
                if part < 3 {
                    w.cross(&Vector3::from(rotated_points[n]))[part]
                } else {
                    w[part - 3]
                }
            }
        });
        assert!(
            (rotated.stiffness * rigid).norm() < 1e-12 * rotated.stiffness.norm() * rigid.norm()
        );
    }
    let mut d = V18::zeros();
    for (n, x) in points.iter().enumerate() {
        d[6 * n] = 0.002 * x[0] + 0.001 * x[1];
        d[6 * n + 1] = -0.003 * x[1] + 0.004 * x[0];
        d[6 * n + 2] = -(0.005 * x[0] * x[0] + 0.007 * x[1] * x[1]) / 2. - 0.002 * x[0] * x[1];
        d[6 * n + 3] = -0.007 * x[1] - 0.002 * x[0];
        d[6 * n + 4] = 0.005 * x[0] + 0.002 * x[1];
        d[6 * n + 5] = (0.004 - 0.001) / 2.;
    }
    let report = kernel.recover("patch".into(), &d);
    for p in &report.samples {
        assert!(
            (Vector3::from(p.membrane_strain) - Vector3::new(0.002, -0.003, 0.005)).norm() < 1e-14
        );
        assert!(
            (Vector3::from(p.curvature_per_m) - Vector3::new(0.005, 0.007, 0.004)).norm() < 1e-14
        );
        assert!(p.drilling_mismatch_rad.abs() < 1e-15);
    }
    let mut rd = V18::zeros();
    for n in 0..3 {
        for part in 0..2 {
            let v = rotation
                * Vector3::new(
                    d[6 * n + 3 * part],
                    d[6 * n + 3 * part + 1],
                    d[6 * n + 3 * part + 2],
                );
            for j in 0..3 {
                rd[6 * n + 3 * part + j] = v[j];
            }
        }
    }
    assert!((d.dot(&(kernel.stiffness * d)) - rd.dot(&(rotated.stiffness * rd))).abs() < 1e-9);
    let energy = report.membrane_energy_kn_m
        + report.bending_energy_kn_m
        + report.drilling_penalty_energy_kn_m;
    assert!((energy - 0.5 * d.dot(&(kernel.stiffness * d))).abs() < 1e-9);
    let eigen = kernel.stiffness.symmetric_eigen().eigenvalues;
    let max = eigen.max();
    assert_eq!(eigen.iter().filter(|v| **v > 1e-10 * max).count(), 12);
    Ok(())
}

#[test]
fn dead_traction_preserves_global_force_and_moment_and_bad_inputs_are_rejected() -> Result<()> {
    let points = [[0.3, -0.8, 2.], [2.1, 1.3, 4.7], [0.9, 2.4, 1.6]];
    let s = section(0.02, 1e-4);
    let kernel = Kernel::new(points, &s)?;
    let traction = Vector3::new(3., -7., 2.);
    let force = kernel.uniform_load(traction.into());
    let mut sum = Vector3::zeros();
    let mut moment = Vector3::zeros();
    for n in 0..3 {
        let f = Vector3::new(force[6 * n], force[6 * n + 1], force[6 * n + 2]);
        sum += f;
        moment += Vector3::from(points[n]).cross(&f);
    }
    let center = points.map(Vector3::from).into_iter().sum::<Vector3<f64>>() / 3.;
    assert!((sum - traction * kernel.area_m2).norm() < 1e-12);
    assert!((moment - center.cross(&(traction * kernel.area_m2))).norm() < 1e-12);
    assert!(Kernel::new([[0.; 3]; 3], &s).is_err());
    assert!(Kernel::new(points, &section(0., 1e-4)).is_err());
    assert!(Kernel::new(points, &section(0.02, 0.)).is_err());
    Ok(())
}

fn strip(divisions: usize, thickness: f64, drilling: f64) -> Result<beams::Configuration> {
    let nodes: Vec<_> = (0..=divisions).flat_map(|i| (0..2).map(move |j|
        json!({"id":format!("n{}",2*i+j), "position_m":[2.*i as f64/divisions as f64,j as f64,0.],"soil_node_set":null}))).collect();
    let plates: Vec<_> = (0..divisions).flat_map(|i| [[2*i,2*i+2,2*i+3],[2*i,2*i+3,2*i+1]].into_iter().enumerate().map(move |(j,cell)|
        json!({"id":format!("p{i}_{j}"),"nodes":cell.map(|n|format!("n{n}")),"section_id":"wall","initially_active":true}))).collect();
    let mut translations = vec![];
    let mut rotations = vec![];
    for n in 0..2 {
        for a in 0..3 {
            let direction: [f64; 3] = std::array::from_fn(|j| if j == a { 1. } else { 0. });
            let c = json!({"node_id":format!("n{n}"),"direction":direction,"value":0.});
            translations.push(c.clone());
            rotations.push(c);
        }
    }
    Ok(serde_json::from_value(
        json!({"schema":"soil-fem-structures/1","nodes":nodes,"sections":[],"members":[],
        "plate_sections":[section(thickness,drilling)],"plates":plates,
        "rotation_length_scale_m":2.,"maximum_condensed_dofs":300,
        "translation_constraints":translations,"rotation_constraints":rotations,
        "nodal_loads":[{"node_id":format!("n{}",2*divisions),"constant_moment_kn_m":[0.,1.,0.]},
                       {"node_id":format!("n{}",2*divisions+1),"constant_moment_kn_m":[0.,1.,0.]}]}),
    )?)
}

#[test]
fn thin_cantilever_plate_matches_closed_form_without_shear_locking_or_drilling_dependence()
-> Result<()> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    for n in [1, 2, 4] {
        for h in [0.02_f64, 0.002] {
            for drilling in [1e-6, 1e-4, 0.01] {
                let c = strip(n, h, drilling)?;
                let model = beams::Model::new(&mesh, c)?;
                let r = model.evaluate(
                    &model.initial_state(),
                    &vec![[0.; 3]; mesh.nodes.len()],
                    &vec![true; mesh.nodes.len()],
                    0.,
                )?;
                let rigidity = 2e8 * h.powi(3) / 12.;
                let curvature = 2. / rigidity;
                for (i, x) in model.positions.iter().enumerate() {
                    assert!(
                        (r.state.node_displacement_m[i][2] + curvature * x[0] * x[0] / 2.).abs()
                            < 1e-7 * curvature
                    );
                    assert!(
                        (r.state.node_rotation_rad[i][1] - curvature * x[0]).abs()
                            < 1e-7 * curvature
                    );
                }
                assert!((r.report.support_moment_origin_kn_m[1] + 2.).abs() < 1e-7);
                assert!(r.report.plates.iter().all(
                    |p| p.drilling_penalty_energy_kn_m < 1e-16 * p.bending_energy_kn_m.max(1.)
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn plate_is_coupled_to_mc_and_mcc_with_installation_cutback_and_checkpoint_history() -> Result<()> {
    for mcc in [false, true] {
        let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
        let base = mesh.nodes.iter().position(|p| *p == [0., 0., 1.]).unwrap();
        mesh.node_sets.insert("PlateBase".into(), vec![base]);
        let top = soil_fem_3d::boundary::surface_owners(&mesh)
            .into_iter()
            .filter(|(face, owners)| {
                owners.len() == 1 && face.iter().all(|&n| mesh.nodes[n][2] == 1.)
            })
            .map(|(face, _)| face)
            .collect();
        mesh.face_sets.insert("top".into(), top);
        let mesh = mesh.elevate_straight()?;
        let ground: Ground = serde_json::from_str(if mcc {
            include_str!("../examples/ground-geostatic-mcc-column.json")
        } else {
            include_str!("../examples/ground-geostatic-mc-column.json")
        })?;
        let mut case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
        let c = json!({"schema":"soil-fem-structures/1","nodes":[
            {"id":"a","position_m":[0.,0.,1.],"soil_node_set":"PlateBase"},
            {"id":"b","position_m":[2.,0.,1.],"soil_node_set":null},
            {"id":"c","position_m":[2.,1.,1.],"soil_node_set":null}],
            "sections":[],"members":[],"plate_sections":[section(0.1,1e-4)],
            "plates":[{"id":"wall","nodes":["a","b","c"],"section_id":"wall","initially_active":false}],
            "rotation_length_scale_m":2.,"maximum_condensed_dofs":30,
            "rotation_constraints":[{"node_id":"a","direction":[1.,0.,0.],"value":0.},{"node_id":"a","direction":[0.,1.,0.],"value":0.},{"node_id":"a","direction":[0.,0.,1.],"value":0.}],
            "nodal_loads":[{"node_id":"b","boundary_force_kn":[0.,0.,-5.]}]});
        case.beam_system = Some(serde_json::from_value(c)?);
        case.settings.initial_increment = 0.25;
        case.stages = serde_json::from_value(json!([
            {"kind":"install_structures","name":"install","members":["wall"],"placement":"cad_new_nodes"},
            {"kind":"boundary","name":"load_plate","target_factor":1.},
            {"kind":"boundary","name":"unload_plate","target_factor":0.},
            {"kind":"remove_structures","name":"remove","members":["wall"]}]))?;
        case.verification_reject_once_stages = vec!["install".into(), "load_plate".into()];
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete, "{:?}", full.failure_reason);
        assert_eq!(full.schema, "soil-fem-stage-result/14");
        for frame in &full.frames {
            assert!(frame.force_balance_error_kn.iter().all(|v| v.abs() < 1e-6));
            assert!(
                frame
                    .global_moment_balance_error_kn_m
                    .unwrap()
                    .iter()
                    .all(|v| v.abs() < 1e-6)
            );
        }
        let loaded = &full.frames[2].beam_results.as_ref().unwrap().plates[0];
        assert!(loaded.bending_energy_kn_m > 0.);
        assert!(full.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_structure_changed
            && a.committed_state_unchanged_on_rejection));
        let initial = full
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .count();
        for count in [initial + 1, initial + 2] {
            let checkpoint = Checkpoint::seal(
                &mesh,
                stages::solve_until(&mesh, ground.clone(), case.clone(), count)?,
            )?;
            let resumed = checkpoint.resume(&mesh, ground.clone(), case.clone(), None)?;
            assert_eq!(
                serde_json::to_value(&full.last_committed)?,
                serde_json::to_value(&resumed.last_committed)?
            );
            assert_eq!(
                serde_json::to_value(&full.frames)?,
                serde_json::to_value(&resumed.frames)?
            );
        }
        let mut invalid = full.last_committed.beams.as_ref().unwrap().clone();
        invalid.plate_reference_global_kinematics.clear();
        let model = beams::Model::new(&mesh, case.beam_system.unwrap())?;
        assert!(model.validate_state(&invalid).is_err());
    }
    Ok(())
}
