//! General initialization plus partial local excavation: never rewind accepted
//! MCC history when a stop arrives immediately before or just after resume.
use soil_fem_3d::{
    Result,
    checkpoint::{self, Checkpoint},
    geostatic::{self, PointSeed, Seed},
    ground::{self, Ground},
    mesh::Mesh,
    stage_control::{Observer, Progress},
    stages::{self, StageCase},
};
use std::cell::Cell;

struct Stop {
    queries: Cell<usize>,
    after: usize,
}
impl Observer for Stop {
    fn on_progress(&self, _: &Progress) -> Result<()> {
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        let n = self.queries.get();
        self.queries.set(n + 1);
        Ok(n >= self.after)
    }
}
fn fixture(projected: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh =
        serde_json::from_str(include_str!("../examples/hardening-column-tet10-mesh.json"))?;
    let cut = mesh
        .tetrahedra
        .iter()
        .position(|t| t.iter().filter(|&&n| mesh.nodes[n][2] == 1.).count() == 3)
        .unwrap();
    mesh.excavation_sets.insert("LocalCut".into(), vec![cut]);
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-mixed-mcc.json"))?;
    for material in &mut ground.materials {
        material.unit_weight_kn_m3 = 1e-8;
        let soil_fem_3d::soil::Model::Mcc(p) = &mut material.model else {
            panic!("MCC required")
        };
        p.critical_state_m = 1.2;
        p.lambda = 0.2;
        p.kappa = 0.04;
        p.poisson = 0.3;
        p.reference_pressure_kpa = 1.;
        p.specific_volume_ncl_at_reference = 3.5;
    }
    for layer in &mut ground.layers {
        layer.k0_x = 1.;
        layer.k0_y = 1.;
        layer.consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 1. };
    }
    ground.excavation_bodies = vec!["LocalCut".into()];
    let points = ground::initialize(&mesh, ground.clone())?
        .point_states
        .into_iter()
        .map(|q| PointSeed {
            element_index: q.element_index,
            local_index: q.local_index,
            material_id: q.material_id,
            position_m: q.centroid_m,
            state: q.constitutive_state.into_state(),
        })
        .collect();
    ground.initialization.method = "equilibrate".into();
    ground.initialization.geostatic = Some(geostatic::Configuration {
        seed: Seed::IntegrationPoints {
            mesh_sha256: checkpoint::hash_json(&mesh)?,
            points,
        },
        initial_boundary_factor: 0.,
        maximum_initial_engineering_strain: 0.05,
        relative_volume_tolerance: 1e-8,
        base_loads: Default::default(),
        hydrostatic_face_sets: vec![],
        verification_reject_once: true,
    });
    ground.layers.clear();
    ground.bodies = mesh
        .regions
        .iter()
        .map(|r| geostatic::Body {
            body: r.body.clone(),
            material_id: r.material_id,
        })
        .collect();
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-mixed-mcc.json"))?;
    case["schema"] = serde_json::json!(if projected {
        "soil-fem-projected-mcc-stages/1"
    } else {
        "soil-fem-mixed-mcc-stages/1"
    });
    case["settings"]["initial_increment"] = serde_json::json!(0.25);
    case["settings"]["maximum_free_dofs"] = serde_json::json!(2000);
    case["settings"]["relative_force_tolerance"] = serde_json::json!(1e-9);
    case["settings"]["relative_linear_tolerance"] = serde_json::json!(1e-11);
    case["boundary"]["constraints"][3]["value_m"] = serde_json::json!(-0.005);
    case["stages"] = serde_json::json!([
        {"kind":"boundary","name":"load","target_factor":1.},
        {"kind":"boundary","name":"unload","target_factor":0.8},
        {"kind":"excavate","name":"cut","bodies":["LocalCut"],"remaining_surface_pressure_kpa":0.},
        {"kind":"boundary","name":"after_cut_unload","target_factor":0.6}
    ]);
    case["pressure_on_original_surface_only"] = serde_json::json!(true);
    case["verification_reject_once_stages"] = serde_json::json!(["cut"]);
    Ok((mesh, ground, serde_json::from_value(case)?))
}
fn scientific_bytes(run: &stages::Run) -> Result<Vec<u8>> {
    fn strip_time(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Object(m) => {
                m.remove("elapsed_seconds");
                for x in m.values_mut() {
                    strip_time(x);
                }
            }
            serde_json::Value::Array(a) => {
                for x in a {
                    strip_time(x);
                }
            }
            _ => {}
        }
    }
    let mut value = serde_json::to_value(run)?;
    strip_time(&mut value);
    // Re-serialization also distinguishes the sign of zero, unlike Value equality.
    Ok(serde_json::to_vec(&value)?)
}
fn prove(projected: bool) -> Result<()> {
    let (mesh, ground, case) = fixture(projected)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    assert!(
        full.attempts
            .iter()
            .any(|a| a.stage == "cut" && !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    let geo = full.geostatic.as_ref().unwrap();
    assert!(geo.complete);
    let count = geo.attempts.iter().filter(|a| a.accepted).count()
        + full
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .position(|a| a.stage == "cut")
            .unwrap()
        + 1;
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
    assert!(
        partial.last_committed.stage_progress > 0. && partial.last_committed.stage_progress < 1.
    );
    assert_eq!(
        partial
            .last_committed
            .active_elements
            .iter()
            .filter(|&&a| !a)
            .count(),
        1
    );
    assert!(
        partial
            .last_committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|x| x.abs() > 1e-7))
    );
    assert!(!partial.excavation_audits.is_empty());
    assert!(
        partial
            .last_committed
            .excavation_release_load_kn
            .iter()
            .any(|x| x.abs() > 1e-8)
    );
    let saved = serde_json::to_vec(&partial)?;
    for after in [0, 1] {
        let clone: stages::Run = serde_json::from_slice(&saved)?;
        let stop = Stop {
            queries: Cell::new(0),
            after,
        };
        let stopped = Checkpoint::seal(&mesh, clone)?.resume_observed(
            &mesh,
            ground.clone(),
            case.clone(),
            &stop,
        )?;
        assert_eq!(
            stopped.failure_reason.as_deref(),
            Some("checkpoint_requested")
        );
        if after == 0 {
            assert!(
                serde_json::to_vec(&stopped)? == saved,
                "early stop altered an accepted byte"
            );
        } else {
            assert_eq!(
                stopped.last_committed.stage_index,
                partial.last_committed.stage_index
            );
            assert!(stopped.last_committed.stage_progress > partial.last_committed.stage_progress);
        }
        assert_eq!(
            stopped.last_committed.active_elements,
            partial.last_committed.active_elements
        );
        let cut = mesh.excavation_sets["LocalCut"][0];
        for (i, q) in mesh.integration_layout()?.into_iter().enumerate() {
            if q.element_index == cut {
                assert!(
                    serde_json::to_vec(&stopped.last_committed.material_states[i])?
                        == serde_json::to_vec(&partial.last_committed.material_states[i])?
                );
            }
        }
        let sealed = serde_json::to_vec(&Checkpoint::seal(&mesh, stopped)?)?;
        let sealed: Checkpoint = serde_json::from_slice(&sealed)?;
        let resumed = sealed.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert!(resumed.complete);
        assert!(
            scientific_bytes(&resumed)? == scientific_bytes(&full)?,
            "complete accepted history changed on restart"
        );
    }
    Ok(())
}
#[test]
fn mixed_mcc_local_cut_early_and_late_stop_restore_every_accepted_history_bit() -> Result<()> {
    prove(false)
}
#[test]
fn projected_mcc_local_cut_early_and_late_stop_restore_every_accepted_history_bit() -> Result<()> {
    prove(true)
}
