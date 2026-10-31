//! H-refinement checks for displacement, stress and a uniform-stress SRM factor.
//!
//! Smooth elasticity goes through `mechanics::solve`. The uniform MC block goes
//! through the accepted-origin epoch search. Neither result is a certificate of
//! general 3D slope convergence or of P6.
use crate::{
    Error, Result,
    boundary::{Conditions, Constraint, NodalLoad, SurfaceLoad},
    elastic, elements,
    epoch_srm_search::{self, Specification as EpochSpec},
    geostatic,
    ground::{Ground, Initialization, Material as GroundMaterial},
    linalg::{GmresOptions, LinearSolver},
    mc, mechanics,
    mesh::{Mesh, Region},
    replay::{PsiPolicy, TensionPolicy},
    soil,
    srm_search::{Outcome, Requirement, Specification as SearchSpec, StopReason},
    stages::{self, Kinematics, Settings, Stage, StageCase},
};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

const SCHEMA: &str = "soil-fem-spatial-convergence/1";
/// Tresca uniaxial strength 2c over the applied 15 kPa traction.
pub const ANALYTICAL_SRM_FACTOR: f64 = 4. / 3.;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Profile {
    /// Nested meshes small enough for the ordinary release test run.
    Ci,
    /// One finer elastic tier and a tighter SRM spacing, for a workstation.
    User,
}

#[derive(Clone, Copy)]
struct Iso {
    young_kpa: f64,
    poisson: f64,
    lambda: f64,
    mu: f64,
    constrained_modulus_kpa: f64,
}

impl Iso {
    fn new(young_kpa: f64, poisson: f64) -> Result<Self> {
        let _ = elastic::constitutive(young_kpa, poisson)?;
        let lambda = young_kpa * poisson / ((1. + poisson) * (1. - 2. * poisson));
        let mu = young_kpa / (2. * (1. + poisson));
        Ok(Self {
            young_kpa,
            poisson,
            lambda,
            mu,
            constrained_modulus_kpa: young_kpa * (1. - poisson)
                / ((1. + poisson) * (1. - 2. * poisson)),
        })
    }
}

#[derive(Serialize)]
pub struct Level {
    pub divisions: usize,
    pub element: String,
    pub h_m: f64,
    pub nodes: usize,
    pub elements: usize,
    pub displacement_l2_m: f64,
    pub stress_l2_kpa: f64,
    pub monitor_m: f64,
    pub monitor_exact_m: f64,
    pub monitor_error_m: f64,
    pub true_relative_free_force_residual: f64,
}

#[derive(Serialize)]
pub struct Series {
    pub id: String,
    pub description: String,
    pub levels: Vec<Level>,
    pub displacement_orders: Vec<Option<f64>>,
    pub stress_orders: Vec<Option<f64>>,
    pub monitor: Extrapolation,
}

#[derive(Serialize)]
pub struct Extrapolation {
    pub quantity: String,
    pub reference: f64,
    pub values: Vec<f64>,
    pub errors: Vec<f64>,
    pub observed_orders: Vec<Option<f64>>,
    pub richardson_estimate: Option<f64>,
    pub richardson_order: Option<f64>,
    pub estimate_minus_reference: Option<f64>,
    pub order_observable: bool,
    pub note: String,
}

#[derive(Serialize)]
pub struct SrmLevel {
    pub divisions: usize,
    pub element: String,
    pub h_m: f64,
    pub nodes: usize,
    pub elements: usize,
    pub attained_equilibrium_factor: Option<f64>,
    pub independently_excluded_factor: Option<f64>,
    pub width: Option<f64>,
    pub midpoint: Option<f64>,
    pub midpoint_minus_analytical: Option<f64>,
    pub stop_reason: String,
    pub requirement: String,
    pub trials: usize,
    pub analytical_factor_inside_bracket: bool,
}

#[derive(Serialize)]
pub struct Report {
    pub schema: String,
    pub profile: String,
    pub p6_complete: bool,
    pub general_3d_spatial_convergence_demonstrated: bool,
    pub smooth_linear_elasticity_orders_demonstrated: bool,
    pub uniform_srm_factor_within_search_spacing: bool,
    pub analytical_uniform_srm_factor: f64,
    pub srm_search_spacing: f64,
    pub elasticity: Vec<Series>,
    pub srm: Vec<SrmLevel>,
    pub srm_extrapolation: Extrapolation,
    pub manufactured_balance_error: f64,
    pub limitations: Vec<String>,
    pub user_command: String,
}

pub fn run(profile: Profile) -> Result<Report> {
    let balance = manufactured_balance_error();
    if balance > 1e-8 {
        return Err(Error::Input(format!(
            "manufactured body force is not in equilibrium with its stress, residual {balance:e}"
        )));
    }
    let (elastic_levels, srm_levels, spacing) = match profile {
        Profile::Ci => (
            vec![
                ("gravity-tet4", false, vec![2, 4, 8]),
                ("gravity-tet10", true, vec![2, 4]),
                ("sine-tet4", false, vec![2, 4, 8]),
                ("sine-tet10", true, vec![2, 4]),
            ],
            vec![(false, 1usize), (false, 2), (true, 1)],
            0.05,
        ),
        Profile::User => (
            vec![
                ("gravity-tet4", false, vec![2, 4, 8, 16]),
                ("gravity-tet10", true, vec![2, 4, 8]),
                ("sine-tet4", false, vec![2, 4, 8, 16]),
                ("sine-tet10", true, vec![2, 4, 8]),
            ],
            vec![
                (false, 1usize),
                (false, 2),
                (false, 4),
                (true, 1),
                (true, 2),
            ],
            0.01,
        ),
    };
    let mut elasticity = Vec::new();
    for (id, tet10, divisions) in elastic_levels {
        elasticity.push(if id.starts_with("gravity") {
            gravity_series(id, tet10, &divisions)?
        } else {
            sine_series(id, tet10, &divisions)?
        });
    }
    let mut srm = Vec::new();
    for (tet10, divisions) in srm_levels {
        srm.push(srm_level(tet10, divisions, spacing)?);
    }
    let smooth = elasticity.iter().all(series_orders_ok);
    let uniform = srm.iter().all(|level| {
        level.analytical_factor_inside_bracket
            && level
                .midpoint_minus_analytical
                .is_some_and(|e| e.abs() <= spacing)
    });
    let srm_extrapolation = extrapolate(
        "uniform MC epoch-SRM bracket midpoint",
        ANALYTICAL_SRM_FACTOR,
        srm.iter()
            .filter(|level| level.element == "tet4")
            .map(|level| level.midpoint.unwrap_or(f64::NAN))
            .collect(),
        spacing,
        "A uniform traction whose linear solution is inside both Tet4 and Tet10 spaces. A flat error at the search spacing is exactness for this patch, not an observed asymptotic order and not a slope Fs.",
    );
    Ok(Report {
        schema: SCHEMA.into(),
        profile: match profile {
            Profile::Ci => "ci",
            Profile::User => "user",
        }
        .into(),
        p6_complete: false,
        general_3d_spatial_convergence_demonstrated: false,
        smooth_linear_elasticity_orders_demonstrated: smooth,
        uniform_srm_factor_within_search_spacing: uniform,
        analytical_uniform_srm_factor: ANALYTICAL_SRM_FACTOR,
        srm_search_spacing: spacing,
        elasticity,
        srm,
        srm_extrapolation,
        manufactured_balance_error: balance,
        limitations: vec![
            "P6 remains incomplete. General 3D spatial convergence of excavation, slopes and non-uniform plastic mechanisms is not demonstrated.".into(),
            "The gravity column is quadratic. Tet10 can represent it, so its error sits at the linear-solver tolerance and does not show an asymptotic order.".into(),
            "The sine body force is applied as consistent nodal loads on the production stiffness. It checks the element spaces; it is not a second constitutive model.".into(),
            "The SRM block is the existing dry φ=0 uniaxial traction with six rigid-body gauges and the epoch search. Its factor 4/3 is uniform-stress exact. Numerical nonconvergence is never used as an upper bound.".into(),
            "The original 2D displacement 5% comparison is not performed or relaxed here.".into(),
        ],
        user_command: "cargo run --release -- spatial-convergence user runs/spatial-convergence-user.json".into(),
    })
}

fn finest_order(orders: &[Option<f64>]) -> Option<f64> {
    orders.iter().rev().find_map(|order| *order)
}

fn series_orders_ok(series: &Series) -> bool {
    let decreasing = series.levels.windows(2).all(|w| {
        w[1].displacement_l2_m < w[0].displacement_l2_m && w[1].stress_l2_kpa < w[0].stress_l2_kpa
    });
    // Quadratic gravity is inside the Tet10 space, so the error is the linear tolerance.
    if series.id.starts_with("gravity") && series.id.ends_with("tet10") {
        return series
            .levels
            .iter()
            .all(|level| level.monitor_error_m < 1e-8 && level.stress_l2_kpa < 1e-6);
    }
    let (disp_band, stress_band) = if series.id.ends_with("tet10") {
        (2.5..3.6, 1.4..2.6)
    } else {
        // A coarse interval can still be pre-asymptotic. The finest interval is the check.
        (1.6..2.4, 0.8..1.3)
    };
    decreasing
        && finest_order(&series.displacement_orders).is_some_and(|order| disp_band.contains(&order))
        && finest_order(&series.stress_orders).is_some_and(|order| stress_band.contains(&order))
        && series.monitor.order_observable
        && series
            .monitor
            .estimate_minus_reference
            .is_some_and(|error| error.abs() < 1e-6)
}

fn gravity_series(id: &str, tet10: bool, divisions: &[usize]) -> Result<Series> {
    let iso = Iso::new(10_000., 0.25)?;
    let gamma = 20.;
    let exact_top = -gamma / (2. * iso.constrained_modulus_kpa);
    let mut levels = Vec::new();
    for &n in divisions {
        let mut mesh = cube(n, tet10)?;
        let supports = rollers(&mesh);
        install(&mut mesh, supports);
        let case = elastic_case(&mesh, iso, [0., 0., -gamma], &[], false)?;
        let solution = mechanics::solve(&mesh, case)?;
        let errors = field_errors(
            &mesh,
            &solution,
            &|x| [0., 0., gravity_uz(x[2], gamma, iso)],
            &|x| gravity_stress(x[2], gamma, iso),
        )?;
        let monitor = top_monitor(&mesh, &solution.displacement_m);
        levels.push(Level {
            divisions: n,
            element: element_name(tet10).into(),
            h_m: 1. / n as f64,
            nodes: mesh.nodes.len(),
            elements: mesh.tetrahedra.len(),
            displacement_l2_m: errors.0,
            stress_l2_kpa: errors.1,
            monitor_m: monitor,
            monitor_exact_m: exact_top,
            monitor_error_m: (monitor - exact_top).abs(),
            true_relative_free_force_residual: solution.true_relative_free_force_residual,
        });
    }
    Ok(finish_series(
        id,
        "Oedometer self-weight. Exact displacement is quadratic and exact stress is linear.",
        levels,
        exact_top,
        "Top-face vertical displacement, shared by every top node.",
    ))
}

fn sine_series(id: &str, tet10: bool, divisions: &[usize]) -> Result<Series> {
    let iso = Iso::new(1_000., 0.3)?;
    let amplitude = 1e-3;
    let mut levels = Vec::new();
    for &n in divisions {
        let mut mesh = cube(n, tet10)?;
        let mut nodal = vec![[0.; 3]; mesh.nodes.len()];
        add_body_force(&mesh, &mut nodal, &|x| sine_body(x, amplitude, iso))?;
        let mut named = BTreeMap::new();
        let loads: Vec<_> = nodal
            .iter()
            .enumerate()
            .filter(|(_, f)| f.iter().any(|v| v.abs() > 0.))
            .map(|(i, f)| {
                let name = format!("load-{i}");
                named.insert(name.clone(), vec![i]);
                NodalLoad {
                    node_set: name,
                    total_force_kn: *f,
                }
            })
            .collect();
        let supports = skin(&mesh, named);
        install(&mut mesh, supports);
        let case = elastic_case(&mesh, iso, [0.; 3], &loads, true)?;
        let solution = mechanics::solve(&mesh, case)?;
        let errors = field_errors(&mesh, &solution, &|x| sine_u(x, amplitude), &|x| {
            sine_stress(x, amplitude, iso)
        })?;
        let monitor = center_monitor(&mesh, &solution.displacement_m);
        let exact = amplitude;
        levels.push(Level {
            divisions: n,
            element: element_name(tet10).into(),
            h_m: 1. / n as f64,
            nodes: mesh.nodes.len(),
            elements: mesh.tetrahedra.len(),
            displacement_l2_m: errors.0,
            stress_l2_kpa: errors.1,
            monitor_m: monitor,
            monitor_exact_m: exact,
            monitor_error_m: (monitor - exact).abs(),
            true_relative_free_force_residual: solution.true_relative_free_force_residual,
        });
    }
    Ok(finish_series(
        id,
        "Unit-cube Dirichlet manufactured solution u_z = A sin(πx) sin(πy) sin(πz), carried by the production stiffness and consistent nodal loads.",
        levels,
        amplitude,
        "Vertical displacement of the node nearest the cube centre.",
    ))
}

fn finish_series(
    id: &str,
    description: &str,
    levels: Vec<Level>,
    reference: f64,
    quantity: &str,
) -> Series {
    let displacement_orders = successive_orders(levels.iter().map(|level| level.displacement_l2_m));
    let stress_orders = successive_orders(levels.iter().map(|level| level.stress_l2_kpa));
    let monitor = extrapolate(
        quantity,
        reference,
        levels.iter().map(|level| level.monitor_m).collect(),
        0.,
        "",
    );
    Series {
        id: id.into(),
        description: description.into(),
        levels,
        displacement_orders,
        stress_orders,
        monitor,
    }
}

fn successive_orders(values: impl Iterator<Item = f64>) -> Vec<Option<f64>> {
    let values: Vec<_> = values.collect();
    values
        .windows(2)
        .map(|w| observed_order(w[0], w[1], 2.))
        .collect()
}

fn observed_order(coarse: f64, fine: f64, ratio: f64) -> Option<f64> {
    if coarse > 0. && fine > 0. && fine < 0.99 * coarse {
        Some((coarse / fine).ln() / ratio.ln())
    } else {
        None
    }
}

fn extrapolate(
    quantity: &str,
    reference: f64,
    values: Vec<f64>,
    floor: f64,
    preset_note: &str,
) -> Extrapolation {
    let errors: Vec<_> = values.iter().map(|v| (v - reference).abs()).collect();
    let observed_orders = successive_orders(errors.iter().copied());
    let floor = floor.max(1e-9);
    let order_observable = observed_orders.iter().any(Option::is_some)
        && errors.last().is_some_and(|error| *error > floor)
        && errors
            .windows(2)
            .any(|pair| pair[0] > floor && pair[1] < 0.7 * pair[0]);
    let (richardson_order, richardson_estimate, estimate_minus_reference, mut note) =
        if order_observable && values.len() >= 2 {
            let order = observed_orders.iter().rev().find_map(|order| *order);
            let mid = values[values.len() - 2];
            let fine = *values.last().unwrap();
            let estimate = order.and_then(|order| richardson(mid, fine, 2., order));
            let delta = estimate.map(|value| value - reference);
            (
                order,
                estimate,
                delta,
                "Observed order from the last error ratio; Richardson uses that order on the two finest values.".into(),
            )
        } else if !errors.is_empty() && errors.iter().all(|error| *error <= floor) {
            let estimate = values.last().copied();
            (
                None,
                estimate,
                estimate.map(|value| value - reference),
                "Successive changes sit at the solver tolerance or the SRM search spacing, so no asymptotic order is identified. The estimate is the finest computed value, and its distance to the reference is the remaining uncertainty.".into(),
            )
        } else {
            (
                None,
                None,
                None,
                "Successive errors do not show a stable positive order. No extrapolated value is reported.".into(),
            )
        };
    if !preset_note.is_empty() {
        note = format!("{preset_note} {note}");
    }
    Extrapolation {
        quantity: quantity.into(),
        reference,
        values,
        errors,
        observed_orders,
        richardson_estimate,
        richardson_order,
        estimate_minus_reference,
        order_observable,
        note,
    }
}

fn richardson(mid: f64, fine: f64, ratio: f64, order: f64) -> Option<f64> {
    let denom = ratio.powf(order) - 1.;
    (denom.is_finite() && denom.abs() > 1e-8).then_some(fine + (fine - mid) / denom)
}

fn srm_level(tet10: bool, divisions: usize, spacing: f64) -> Result<SrmLevel> {
    let mesh = cube(divisions, tet10)?;
    let (ground, case) = srm_inputs(&mesh)?;
    let origin = stages::solve(&mesh, ground, case)?;
    if !origin.complete {
        return Err(Error::Solve(
            "unreduced uniform SRM block did not finish; no factor is accepted".into(),
        ));
    }
    let spec = EpochSpec {
        schema: "soil-fem-epoch-srm-search/1".into(),
        search: SearchSpec {
            schema: "soil-fem-srm-search/1".into(),
            material_ids: vec![1],
            psi_policy: PsiPolicy::Zero,
            tension_policy: TensionPolicy::Keep,
            maximum_factor: 2.,
            required_factor: 1.1,
            absolute_factor_tolerance: spacing,
            maximum_trials: if spacing <= 0.01 { 16 } else { 10 },
            through_stages: 1,
            virtual_displacements_m: vec![],
        },
        initial_increment: 0.25,
        minimum_increment: 1. / 4096.,
        maximum_attempts: 128,
    };
    let report = epoch_srm_search::run(&mesh, origin, spec)?;
    let attained = report
        .interval
        .as_ref()
        .map(|b| b.attained_equilibrium_factor);
    let excluded = report
        .interval
        .as_ref()
        .map(|b| b.independently_excluded_factor);
    let width = report.interval.as_ref().map(|b| b.width);
    let midpoint = match (attained, excluded) {
        (Some(lo), Some(hi)) => Some(0.5 * (lo + hi)),
        _ => None,
    };
    let inside = match (attained, excluded) {
        (Some(lo), Some(hi)) => lo < ANALYTICAL_SRM_FACTOR && ANALYTICAL_SRM_FACTOR < hi,
        _ => false,
    };
    if report
        .trials
        .iter()
        .any(|trial| trial.evidence.outcome == Outcome::Indeterminate)
        && report.interval.is_none()
    {
        return Err(Error::Solve(
            "uniform SRM search stopped indeterminate before an independent bracket; this is recorded as a solve failure, not as an upper factor".into(),
        ));
    }
    Ok(SrmLevel {
        divisions,
        element: element_name(tet10).into(),
        h_m: 1. / divisions as f64,
        nodes: mesh.nodes.len(),
        elements: mesh.tetrahedra.len(),
        attained_equilibrium_factor: attained,
        independently_excluded_factor: excluded,
        width,
        midpoint,
        midpoint_minus_analytical: midpoint.map(|q| q - ANALYTICAL_SRM_FACTOR),
        stop_reason: stop_name(report.stop_reason).into(),
        requirement: requirement_name(report.requirement).into(),
        trials: report.trials.len(),
        analytical_factor_inside_bracket: inside,
    })
}

fn srm_inputs(mesh: &Mesh) -> Result<(Ground, StageCase)> {
    let material = mc::Parameters {
        young_kpa: 100_000.,
        poisson: 0.3,
        cohesion_kpa: 10.,
        friction_deg: 0.,
        dilation_deg: 0.,
        tension_cutoff_kpa: None,
        flow: mc::Flow::NonAssociated,
    };
    let zero = soil::State::MohrCoulomb(mc::State {
        effective_stress_compression_positive_kpa: [0.; 6],
        plastic_strain_compression_positive: [0.; 6],
        void_ratio: 0.7,
        equivalent_plastic_strain: 0.,
        plastic_work_kpa: 0.,
        last_return: mc::Return::Elastic,
    });
    let mut states = BTreeMap::new();
    states.insert("Cube".into(), zero);
    let ground = Ground {
        schema: "soil-fem-ground/1".into(),
        units: "m-kN-kPa".into(),
        initialization: Initialization {
            method: "equilibrate".into(),
            surface_surcharge_kpa: 0.,
            water_level_m: None,
            water_unit_weight_kn_m3: 9.80665,
            relative_equilibrium_tolerance: 1e-9,
            geostatic: Some(geostatic::Configuration {
                seed: geostatic::Seed::BodyStates { states },
                initial_boundary_factor: 0.,
                maximum_initial_engineering_strain: 0.05,
                relative_volume_tolerance: 1e-8,
                base_loads: Default::default(),
                hydrostatic_face_sets: vec![],
                verification_reject_once: false,
            }),
        },
        materials: vec![GroundMaterial {
            id: 1,
            name: "uniform spatial SRM".into(),
            unit_weight_kn_m3: 1e-10,
            model: soil::Model::MohrCoulomb(material),
        }],
        construction_material_ids: vec![],
        layers: vec![],
        bodies: vec![geostatic::Body {
            body: "Cube".into(),
            material_id: 1,
        }],
        excavation_bodies: vec![],
    };
    let case = StageCase {
        joint_system: None,
        anchor_system: None,
        beam_system: None,
        schema: "soil-fem-stages/1".into(),
        settings: Settings {
            kinematics: Kinematics::Explicit,
            relative_force_tolerance: 1e-8,
            relative_linear_tolerance: 1e-10,
            maximum_newton_iterations: 35,
            initial_increment: 1.,
            minimum_increment: 1. / 512.,
            material_maximum_substep: 1e-3,
            mcc_integration: None,
            tangent_perturbation: 1e-7,
            maximum_free_dofs: 100_000,
            linear_solver: sparse(),
        },
        boundary: Some(srm_boundary(mesh)),
        monitors: BTreeMap::new(),
        monitor_sets: BTreeMap::new(),
        pressure_on_original_surface_only: true,
        stages: vec![Stage::Boundary {
            name: "load".into(),
            target_factor: 1.,
            initial_increment: Some(1.),
        }],
        initially_inactive_bodies: vec![],
        verification_reject_once_stages: vec![],
    };
    Ok((ground, case))
}

fn stop_name(reason: StopReason) -> &'static str {
    match reason {
        StopReason::FactorTolerance => "factor_tolerance",
        StopReason::NoPhysicalUpper => "no_physical_upper",
        StopReason::IndeterminateTrial => "indeterminate_trial",
        StopReason::TrialLimit => "trial_limit",
        StopReason::ContradictoryEvidence => "contradictory_evidence",
    }
}

fn requirement_name(requirement: Requirement) -> &'static str {
    match requirement {
        Requirement::EquilibriumRequirementMet => "equilibrium_requirement_met",
        Requirement::IndependentlyExcluded => "independently_excluded",
        Requirement::Indeterminate => "indeterminate",
    }
}

fn srm_boundary(_mesh: &Mesh) -> Conditions {
    let gauge = |set: &str, direction: [f64; 3]| Constraint {
        node_set: set.into(),
        direction,
        value_m: 0.,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    };
    let traction = |set: &str, z: f64| SurfaceLoad::Traction {
        face_set: set.into(),
        traction_kpa: [0., 0., z],
        gradient_kpa_per_m: [[0.; 3]; 3],
        reference_point_m: [0.; 3],
    };
    Conditions {
        constraints: vec![
            gauge("origin", [1., 0., 0.]),
            gauge("origin", [0., 1., 0.]),
            gauge("origin", [0., 0., 1.]),
            gauge("x-gauge", [0., 1., 0.]),
            gauge("x-gauge", [0., 0., 1.]),
            gauge("y-gauge", [0., 0., 1.]),
        ],
        nodal_loads: vec![],
        surface_loads: vec![traction("bottom", 15.), traction("top", -15.)],
    }
}

fn install(mesh: &mut Mesh, sets: BTreeMap<String, Vec<usize>>) {
    for (name, nodes) in sets {
        mesh.node_sets.insert(name, nodes);
    }
}

fn elastic_case(
    _mesh: &Mesh,
    iso: Iso,
    body: [f64; 3],
    loads: &[NodalLoad],
    fix_all_components: bool,
) -> Result<mechanics::Case> {
    let constraints = if fix_all_components {
        ["x", "y", "z"]
            .into_iter()
            .enumerate()
            .map(|(axis, _)| {
                let mut direction = [0.; 3];
                direction[axis] = 1.;
                Constraint {
                    node_set: "skin".into(),
                    direction,
                    value_m: 0.,
                    gradient: [0.; 3],
                    reference_point_m: [0.; 3],
                }
            })
            .collect()
    } else {
        vec![
            axis_constraint("bottom", 2),
            axis_constraint("xmin", 0),
            axis_constraint("xmax", 0),
            axis_constraint("ymin", 1),
            axis_constraint("ymax", 1),
        ]
    };
    Ok(mechanics::Case {
        schema: "soil-fem-mechanics/1".into(),
        units: "m-kN-kPa".into(),
        materials: vec![mechanics::Material {
            id: 1,
            young_kpa: iso.young_kpa,
            poisson: iso.poisson,
            body_force_kn_m3: body,
        }],
        boundary: Conditions {
            constraints,
            nodal_loads: loads.to_vec(),
            surface_loads: vec![],
        },
        relative_linear_tolerance: 1e-10,
        maximum_free_dofs: 100_000,
        linear_solver: sparse(),
    })
}

fn axis_constraint(set: &str, axis: usize) -> Constraint {
    let mut direction = [0.; 3];
    direction[axis] = 1.;
    Constraint {
        node_set: set.into(),
        direction,
        value_m: 0.,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    }
}

fn sparse() -> LinearSolver {
    LinearSolver::SparseGmresIlu0 {
        options: GmresOptions::default(),
    }
}

fn rollers(mesh: &Mesh) -> BTreeMap<String, Vec<usize>> {
    let mut sets = BTreeMap::new();
    let mut bottom = Vec::new();
    let mut xmin = Vec::new();
    let mut xmax = Vec::new();
    let mut ymin = Vec::new();
    let mut ymax = Vec::new();
    for (i, p) in mesh.nodes.iter().enumerate() {
        if p[2].abs() <= 1e-12 {
            bottom.push(i);
        }
        if p[0].abs() <= 1e-12 {
            xmin.push(i);
        }
        if (p[0] - 1.).abs() <= 1e-12 {
            xmax.push(i);
        }
        if p[1].abs() <= 1e-12 {
            ymin.push(i);
        }
        if (p[1] - 1.).abs() <= 1e-12 {
            ymax.push(i);
        }
    }
    sets.insert("bottom".into(), bottom);
    sets.insert("xmin".into(), xmin);
    sets.insert("xmax".into(), xmax);
    sets.insert("ymin".into(), ymin);
    sets.insert("ymax".into(), ymax);
    sets
}

fn skin(mesh: &Mesh, mut sets: BTreeMap<String, Vec<usize>>) -> BTreeMap<String, Vec<usize>> {
    let skin = mesh
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            (p.iter()
                .any(|v| v.abs() <= 1e-12 || (v - 1.).abs() <= 1e-12))
            .then_some(i)
        })
        .collect();
    sets.insert("skin".into(), skin);
    sets
}

fn cube(divisions: usize, tet10: bool) -> Result<Mesh> {
    if !(1..=32).contains(&divisions) {
        return Err(Error::Input("division count must be 1..=32".into()));
    }
    let n = divisions;
    let mut nodes = Vec::new();
    for iz in 0..=n {
        for iy in 0..=n {
            for ix in 0..=n {
                nodes.push([
                    ix as f64 / n as f64,
                    iy as f64 / n as f64,
                    iz as f64 / n as f64,
                ]);
            }
        }
    }
    let id = |ix: usize, iy: usize, iz: usize| ix + (n + 1) * (iy + (n + 1) * iz);
    let mut tetrahedra = Vec::new();
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let c = [
                    id(ix, iy, iz),
                    id(ix + 1, iy, iz),
                    id(ix + 1, iy + 1, iz),
                    id(ix, iy + 1, iz),
                    id(ix, iy, iz + 1),
                    id(ix + 1, iy, iz + 1),
                    id(ix + 1, iy + 1, iz + 1),
                    id(ix, iy + 1, iz + 1),
                ];
                for mut tet in [
                    [c[0], c[1], c[2], c[6]],
                    [c[0], c[2], c[3], c[6]],
                    [c[0], c[3], c[7], c[6]],
                    [c[0], c[7], c[4], c[6]],
                    [c[0], c[4], c[5], c[6]],
                    [c[0], c[5], c[1], c[6]],
                ] {
                    let p = tet.map(|i| nodes[i]);
                    if crate::mesh::geometry(&p)?.0 <= 0. {
                        tet.swap(1, 2);
                    }
                    tetrahedra.push(tet);
                }
            }
        }
    }
    let (bottom, top) = horizontal_faces(&nodes, &tetrahedra);
    let mut node_sets = BTreeMap::new();
    node_sets.insert("origin".into(), vec![node_index(&nodes, [0., 0., 0.])]);
    node_sets.insert("x-gauge".into(), vec![node_index(&nodes, [1., 0., 0.])]);
    node_sets.insert("y-gauge".into(), vec![node_index(&nodes, [0., 1., 0.])]);
    let mut face_sets = BTreeMap::new();
    face_sets.insert("bottom".into(), bottom);
    face_sets.insert("top".into(), top);
    let hash = format!("{:x}", Sha256::digest(b"spatial-convergence-cube"));
    let mesh = Mesh {
        schema: "soil-fem-mesh/2".into(),
        units: "m-kN-kPa".into(),
        axis: "z_up".into(),
        crs: "LOCAL".into(),
        origin_m: [0.; 3],
        nodes,
        tetrahedra,
        tetrahedra10: vec![],
        quadrature_order: None,
        material_ids: vec![1; n * n * n * 6],
        node_sets,
        face_sets,
        regions: vec![Region {
            id: 1,
            body: "Cube".into(),
            material_id: 1,
            bounds_m: [[0.; 3], [1.; 3]],
            volume_m3: 1.,
            source_sha256: hash,
        }],
        region_ids: vec![1; n * n * n * 6],
        excavation_sets: BTreeMap::new(),
    };
    let volume = mesh.validate()?;
    if (volume - 1.).abs() > 1e-9 {
        return Err(Error::Input(format!("cube volume {volume} is not 1")));
    }
    if !tet10 {
        return Ok(mesh);
    }
    let mut elevated = mesh.elevate_straight()?;
    elevated.quadrature_order = Some(5);
    elevated.validate()?;
    Ok(elevated)
}

fn horizontal_faces(
    nodes: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
) -> (Vec<[usize; 3]>, Vec<[usize; 3]>) {
    let mut owners: BTreeMap<[usize; 3], usize> = BTreeMap::new();
    for tet in tetrahedra {
        for local in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut face = [tet[local[0]], tet[local[1]], tet[local[2]]];
            face.sort_unstable();
            *owners.entry(face).or_default() += 1;
        }
    }
    let mut bottom = Vec::new();
    let mut top = Vec::new();
    for (face, count) in owners {
        if count != 1 {
            continue;
        }
        let z: Vec<_> = face.iter().map(|&i| nodes[i][2]).collect();
        if z.iter().all(|z| z.abs() <= 1e-12) {
            bottom.push(face);
        } else if z.iter().all(|z| (z - 1.).abs() <= 1e-12) {
            top.push(face);
        }
    }
    (bottom, top)
}

fn node_index(nodes: &[[f64; 3]], target: [f64; 3]) -> usize {
    nodes
        .iter()
        .position(|p| p.iter().zip(target).all(|(a, b)| (a - b).abs() <= 1e-12))
        .expect("gauge node")
}

fn element_name(tet10: bool) -> &'static str {
    if tet10 { "tet10" } else { "tet4" }
}

fn add_body_force(
    mesh: &Mesh,
    force: &mut [[f64; 3]],
    body: &impl Fn([f64; 3]) -> [f64; 3],
) -> Result<()> {
    for ei in 0..mesh.tetrahedra.len() {
        let nodes = mesh.cell_nodes(ei);
        if nodes.len() == 4 {
            let corners = mesh.points(&mesh.tetrahedra[ei]);
            for &(r, w) in elements::volume_rule(5)? {
                let (x, det, shape) = map_tet4(&corners, r)?;
                let b = body(x);
                let weight = w * det;
                for i in 0..4 {
                    for a in 0..3 {
                        force[nodes[i]][a] += weight * shape[i] * b[a];
                    }
                }
            }
        } else {
            for q in mesh.integration(ei)? {
                let b = body(q.position);
                for (i, &node) in nodes.iter().enumerate() {
                    for a in 0..3 {
                        force[node][a] += q.weight_m3 * q.shape[i] * b[a];
                    }
                }
            }
        }
    }
    Ok(())
}

fn map_tet4(points: &[[f64; 3]; 4], r: [f64; 3]) -> Result<([f64; 3], f64, [f64; 4])> {
    let shape = [1. - r[0] - r[1] - r[2], r[0], r[1], r[2]];
    let mut x = [0.; 3];
    for i in 0..4 {
        for a in 0..3 {
            x[a] += shape[i] * points[i][a];
        }
    }
    let (volume, _) = crate::mesh::geometry(points)?;
    let det = 6. * volume;
    Ok((x, det, shape))
}

fn field_errors(
    mesh: &Mesh,
    solution: &mechanics::Solution,
    displacement: &impl Fn([f64; 3]) -> [f64; 3],
    stress: &impl Fn([f64; 3]) -> [f64; 6],
) -> Result<(f64, f64)> {
    let mut disp = 0.;
    let mut sig = 0.;
    let mut volume = 0.;
    if mesh.tetrahedra10.is_empty() {
        let iso = &solution.case.materials[0];
        let d = elastic::constitutive(iso.young_kpa, iso.poisson)?;
        for tet in &mesh.tetrahedra {
            let corners = mesh.points(tet);
            let (_, b, _) = elastic::element(&corners, &d)?;
            let ue = std::array::from_fn(|k| solution.displacement_m[tet[k / 3]][k % 3]);
            let fe = elastic::stress(&b, &d, &ue);
            for &(r, w) in elements::volume_rule(4)? {
                let (x, det, shape) = map_tet4(&corners, r)?;
                let weight = w * det;
                let u: [f64; 3] = std::array::from_fn(|a| {
                    (0..4)
                        .map(|i| shape[i] * solution.displacement_m[tet[i]][a])
                        .sum()
                });
                let exact_u = displacement(x);
                let exact_s = stress(x);
                disp += weight * (0..3).map(|a| (u[a] - exact_u[a]).powi(2)).sum::<f64>();
                sig += weight * (0..6).map(|a| (fe[a] - exact_s[a]).powi(2)).sum::<f64>();
                volume += weight;
            }
        }
    } else {
        for (point, fe) in solution
            .integration_points
            .iter()
            .zip(&solution.stress_tension_positive_kpa)
        {
            let nodes = mesh.cell_nodes(point.element_index);
            let u: [f64; 3] = std::array::from_fn(|a| {
                nodes
                    .iter()
                    .zip(&point.shape_functions)
                    .map(|(&n, s)| s * solution.displacement_m[n][a])
                    .sum()
            });
            let exact_u = displacement(point.position_m);
            let exact_s = stress(point.position_m);
            let weight = point.weight_m3;
            disp += weight * (0..3).map(|a| (u[a] - exact_u[a]).powi(2)).sum::<f64>();
            sig += weight * (0..6).map(|a| (fe[a] - exact_s[a]).powi(2)).sum::<f64>();
            volume += weight;
        }
    }
    if volume <= 0. {
        return Err(Error::Input("error norm has no positive volume".into()));
    }
    Ok(((disp / volume).sqrt(), (sig / volume).sqrt()))
}

fn top_monitor(mesh: &Mesh, displacement: &[[f64; 3]]) -> f64 {
    displacement
        .iter()
        .zip(&mesh.nodes)
        .filter(|(_, p)| (p[2] - 1.).abs() <= 1e-12)
        .map(|(u, _)| u[2])
        .sum::<f64>()
        / mesh
            .nodes
            .iter()
            .filter(|p| (p[2] - 1.).abs() <= 1e-12)
            .count() as f64
}

fn center_monitor(mesh: &Mesh, displacement: &[[f64; 3]]) -> f64 {
    let (index, _) = mesh
        .nodes
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let da = a.iter().map(|v| (v - 0.5).powi(2)).sum::<f64>();
            let db = b.iter().map(|v| (v - 0.5).powi(2)).sum::<f64>();
            da.total_cmp(&db)
        })
        .expect("node");
    displacement[index][2]
}

fn gravity_uz(z: f64, gamma: f64, iso: Iso) -> f64 {
    -gamma * (z - 0.5 * z * z) / iso.constrained_modulus_kpa
}

fn gravity_stress(z: f64, gamma: f64, iso: Iso) -> [f64; 6] {
    let vertical = -gamma * (1. - z);
    let horizontal = iso.poisson / (1. - iso.poisson) * vertical;
    [horizontal, horizontal, vertical, 0., 0., 0.]
}

fn sine_u(x: [f64; 3], amplitude: f64) -> [f64; 3] {
    [
        0.,
        0.,
        amplitude
            * (x[0] * std::f64::consts::PI).sin()
            * (x[1] * std::f64::consts::PI).sin()
            * (x[2] * std::f64::consts::PI).sin(),
    ]
}

fn sine_stress(x: [f64; 3], amplitude: f64, iso: Iso) -> [f64; 6] {
    let p = std::f64::consts::PI;
    let sx = (p * x[0]).sin();
    let sy = (p * x[1]).sin();
    let sz = (p * x[2]).sin();
    let cx = (p * x[0]).cos();
    let cy = (p * x[1]).cos();
    let cz = (p * x[2]).cos();
    let trace = amplitude * p * sx * sy * cz;
    let gamma_yz = amplitude * p * sx * cy * sz;
    let gamma_zx = amplitude * p * cx * sy * sz;
    [
        iso.lambda * trace,
        iso.lambda * trace,
        (iso.lambda + 2. * iso.mu) * trace,
        0.,
        iso.mu * gamma_yz,
        iso.mu * gamma_zx,
    ]
}

fn sine_body(x: [f64; 3], amplitude: f64, iso: Iso) -> [f64; 3] {
    let p = std::f64::consts::PI;
    let sx = (p * x[0]).sin();
    let sy = (p * x[1]).sin();
    let sz = (p * x[2]).sin();
    let cx = (p * x[0]).cos();
    let cy = (p * x[1]).cos();
    let cz = (p * x[2]).cos();
    let scale = amplitude * p * p;
    [
        -(iso.lambda + iso.mu) * scale * cx * sy * cz,
        -(iso.lambda + iso.mu) * scale * sx * cy * cz,
        (iso.lambda + 4. * iso.mu) * scale * sx * sy * sz,
    ]
}

pub fn manufactured_balance_error() -> f64 {
    let iso = Iso::new(1_000., 0.3).expect("material");
    let amplitude = 1e-3;
    let mut worst = 0_f64;
    for x in [[0.2, 0.3, 0.4], [0.5, 0.5, 0.5], [0.7, 0.2, 0.8]] {
        let h = 1e-6;
        let mut numerical = [0.; 3];
        for a in 0..3 {
            let mut plus = x;
            let mut minus = x;
            plus[a] += h;
            minus[a] -= h;
            let sp = sine_stress(plus, amplitude, iso);
            let sm = sine_stress(minus, amplitude, iso);
            // σxx,x + σxy,y + σzx,z and cyclic. Voigt is xx,yy,zz,xy,yz,zx.
            let columns = [[0, 3, 5], [3, 1, 4], [5, 4, 2]];
            for (row, cols) in columns.iter().enumerate() {
                numerical[row] += (sp[cols[a]] - sm[cols[a]]) / (2. * h);
            }
        }
        let body = sine_body(x, amplitude, iso);
        let scale = body.iter().map(|v| v.abs()).fold(1., f64::max);
        for a in 0..3 {
            worst = worst.max((numerical[a] + body[a]).abs() / scale);
        }
    }
    worst
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::{self, Case};

    #[test]
    fn manufactured_stress_matches_its_body_force() {
        assert!(manufactured_balance_error() < 1e-8);
    }

    #[test]
    fn existing_single_tet_settlement_is_unchanged() -> Result<()> {
        let mesh = Mesh {
            schema: "soil-fem-mesh/1".into(),
            units: "m-kN-kPa".into(),
            axis: "z_up".into(),
            crs: "LOCAL".into(),
            origin_m: [0.; 3],
            nodes: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            tetrahedra: vec![[0, 1, 2, 3]],
            tetrahedra10: vec![],
            quadrature_order: None,
            material_ids: vec![1],
            node_sets: BTreeMap::from([("base".into(), vec![0, 1, 2])]),
            face_sets: BTreeMap::new(),
            regions: vec![],
            region_ids: vec![],
            excavation_sets: BTreeMap::new(),
        };
        let case = Case {
            schema: "soil-fem-elastic-case/1".into(),
            young_kpa: 1000.,
            poisson: 0.25,
            unit_weight_kn_m3: 20.,
            fixed_node_set: "base".into(),
            relative_tolerance: 1e-10,
            max_iterations: 100,
        };
        let solution = solver::solve(&mesh, &case)?;
        assert!((solution.displacement_m[3][2] + 1. / 240.).abs() < 1e-12);
        assert!((solution.support_reaction_kn[2] - 20. / 6.).abs() < 1e-12);
        assert!(
            solution
                .equilibrium_error_kn
                .iter()
                .all(|x| x.abs() < 1e-12)
        );
        Ok(())
    }

    #[test]
    fn ci_profile_reports_orders_without_closing_p6() -> Result<()> {
        let first = gravity_series("gravity-tet4", false, &[2])?;
        let second = gravity_series("gravity-tet4", false, &[2])?;
        assert_eq!(
            first.levels[0].monitor_m.to_bits(),
            second.levels[0].monitor_m.to_bits()
        );
        let report = run(Profile::Ci)?;
        let summary = serde_json::to_string_pretty(&report)?;
        assert!(!report.p6_complete);
        assert!(!report.general_3d_spatial_convergence_demonstrated);
        assert!(
            report.smooth_linear_elasticity_orders_demonstrated,
            "{summary}"
        );
        assert!(report.uniform_srm_factor_within_search_spacing, "{summary}");
        assert!(report.manufactured_balance_error < 1e-8);
        let gravity = report
            .elasticity
            .iter()
            .find(|s| s.id == "gravity-tet4")
            .unwrap();
        let sine = report
            .elasticity
            .iter()
            .find(|s| s.id == "sine-tet4")
            .unwrap();
        assert!((1.6..2.4).contains(&finest_order(&gravity.displacement_orders).unwrap()));
        assert!((0.8..1.3).contains(&finest_order(&gravity.stress_orders).unwrap()));
        assert!((1.6..2.4).contains(&finest_order(&sine.displacement_orders).unwrap()));
        assert!((0.8..1.3).contains(&finest_order(&sine.stress_orders).unwrap()));
        let sine10 = report
            .elasticity
            .iter()
            .find(|series| series.id == "sine-tet10")
            .unwrap();
        assert!(finest_order(&sine10.displacement_orders).unwrap() > 2.5);
        assert!(finest_order(&sine10.stress_orders).unwrap() > 1.4);
        let gravity10 = report
            .elasticity
            .iter()
            .find(|series| series.id == "gravity-tet10")
            .unwrap();
        assert!(
            gravity10
                .levels
                .iter()
                .all(|level| level.stress_l2_kpa < 1e-6)
        );
        assert!(!gravity10.monitor.order_observable);
        for level in &report.srm {
            assert!(level.analytical_factor_inside_bracket);
            assert!(
                level
                    .width
                    .is_some_and(|width| width <= report.srm_search_spacing + 1e-12)
            );
            assert!(level.midpoint_minus_analytical.unwrap().abs() <= report.srm_search_spacing);
        }
        let midpoints: Vec<_> = report
            .srm
            .iter()
            .map(|level| level.midpoint.unwrap())
            .collect();
        assert!(
            midpoints
                .windows(2)
                .all(|pair| (pair[0] - pair[1]).abs() <= 1e-12)
        );
        for series in [gravity, sine, sine10] {
            let error = series.monitor.estimate_minus_reference.unwrap();
            assert!(
                error.abs() < 1e-6,
                "{} richardson error {error:e}",
                series.id
            );
        }
        Ok(())
    }
}
