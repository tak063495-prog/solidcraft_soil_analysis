//! Plane-strain 3D displacement check against the original 2D column.
//!
//! Depth is fixed. Nested cross-section meshes go through `stages::solve`.
//! The crest monitor uses the same Richardson procedure as
//! `spatial_convergence::extrapolate`. The 5% limit is the archived slope
//! protocol value and is not relaxed. A pass on this column does not close
//! the Griffiths slope comparison or P6.

use crate::{
    Error, Result,
    boundary::{Conditions, Constraint},
    geostatic,
    ground::{Ground, Initialization, Material as GroundMaterial},
    linalg::{GmresOptions, LinearSolver},
    mc,
    mesh::{Mesh, Region},
    soil, spatial_convergence,
    spatial_convergence::Extrapolation,
    stages::{self, Kinematics, Settings, Stage, StageCase},
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

const SCHEMA: &str = "soil-fem-plane-strain-2d/1";
const SERIES_SCHEMA: &str = "soil-fem-plane-strain-2d-series/1";
/// Archived predeclared gate. `displacement_sample_relative_L2` in the
/// Griffiths slope protocol. Not a new tolerance.
pub const DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT: f64 = 0.05;
/// Companion gate from the same protocol. Not used in place of the 5% gate.
pub const MAXIMUM_DISPLACEMENT_ERROR_OVER_REFERENCE_MAX: f64 = 0.10;
/// Archived `depth_displacement_variation_over_max`. Not relaxed, and not a
/// substitute for the displacement L2 gate.
pub const DEPTH_DISPLACEMENT_VARIATION_LIMIT: f64 = 0.01;
/// `docs/2D_COMPARISON.ja.md`: maximum nodal error of the original Q8 run
/// against the closed form. The printed crest values stay inside this band.
pub const ORIGINAL_2D_THEORY_DISPLACEMENT_ERROR_M: f64 = 4.14e-13;

/// Printed original 2D crest displacements, metres.
/// Gravity, loaded, unloaded. From `docs/2D_COMPARISON.ja.md`.
pub const ORIGINAL_2D_CREST_M: [f64; 3] = [-0.333333333e-3, -2.700167230e-3, -0.700167230e-3];

const STAGE_NAMES: [&str; 3] = ["gravity", "loaded", "unloaded"];

/// Middle-depth relative L2 from the archived failure files. These runs are
/// not an in-plane h series, so they are not extrapolated.
pub const ARCHIVED_SLOPE_MIDDLE_RELATIVE_L2: [(&str, f64); 4] = [
    ("xy1-ny1", 0.1558199870691725),
    ("xy1-ny2", 0.15885691929516224),
    ("xy1-ny1-global64", 0.15593219845164077),
    ("xy1-ny1-local-half", 0.15672391195315055),
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Verdict {
    Pass,
    Fail,
    Indeterminate,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Profile {
    Ci,
    User,
}

#[derive(Clone, Copy)]
struct Column {
    width_m: f64,
    depth_m: f64,
    height_m: f64,
    young_kpa: f64,
    poisson: f64,
    unit_weight_kn_m3: f64,
    cohesion_kpa: f64,
    friction_deg: f64,
    pressure_kpa: f64,
}

fn original_column() -> Column {
    Column {
        width_m: 0.2,
        depth_m: 1.0,
        height_m: 0.2,
        young_kpa: 1000.,
        poisson: 0.25,
        unit_weight_kn_m3: 20.,
        cohesion_kpa: 2.,
        friction_deg: 10.,
        pressure_kpa: 12.,
    }
}

/// Closed form matched by the original 2D Q8 column. Compression-positive
/// internals; the returned displacement is z-up.
pub fn original_2d_vertical_displacement_m(z: f64, stage: usize) -> Result<f64> {
    let column = original_column();
    if !(0.0..=column.height_m).contains(&z) || stage > 2 {
        return Err(Error::Input(
            "column displacement sample is outside the original 2D column".into(),
        ));
    }
    let sine = column.friction_deg.to_radians().sin();
    let cosine = column.friction_deg.to_radians().cos();
    let bulk = column.young_kpa / (3. * (1. - 2. * column.poisson));
    let constrained = column.young_kpa * (1. - column.poisson)
        / ((1. + column.poisson) * (1. - 2. * column.poisson));
    let edge = 3. * (1. + sine) / (3. - sine);
    let intercept = 4. * column.cohesion_kpa * cosine / (3. - sine);
    let height = column.height_m;
    let gravity = -column.unit_weight_kn_m3 * (height * z - 0.5 * z * z) / constrained;
    let loaded = -((column.pressure_kpa + column.unit_weight_kn_m3 * height - intercept) * z
        - 0.5 * column.unit_weight_kn_m3 * z * z)
        / (edge * bulk);
    let unloaded = loaded + column.pressure_kpa * z / constrained;
    Ok([gravity, loaded, unloaded][stage])
}

pub fn documented_crest_agrees_with_closed_form() -> Result<bool> {
    let height = original_column().height_m;
    for (stage, printed) in ORIGINAL_2D_CREST_M.iter().enumerate() {
        let closed = original_2d_vertical_displacement_m(height, stage)?;
        if (closed - printed).abs() > ORIGINAL_2D_THEORY_DISPLACEMENT_ERROR_M {
            return Ok(false);
        }
    }
    Ok(true)
}

#[derive(Serialize)]
pub struct FieldSample {
    pub displacement_sample_relative_l2: f64,
    pub maximum_displacement_error_over_reference_max: f64,
    pub depth_displacement_variation_over_reference_max: f64,
    pub maximum_abs_depth_displacement_m: f64,
    pub maximum_abs_lateral_displacement_m: f64,
    pub crest_mean_uz_m: f64,
    pub crest_reference_m: f64,
    pub relative_l2_verdict: Verdict,
    pub maximum_ratio_verdict: Verdict,
    pub depth_variation_verdict: Verdict,
}

#[derive(Serialize)]
pub struct Level {
    pub divisions: usize,
    pub element: String,
    pub h_m: f64,
    pub nodes: usize,
    pub elements: usize,
    pub stages: Vec<FieldSample>,
}

#[derive(Serialize)]
pub struct ArchivedSlope {
    pub case: String,
    pub middle_depth_displacement_sample_relative_l2: f64,
    pub five_percent_verdict: Verdict,
    pub in_plane_h_series: bool,
    pub continuum_five_percent_verdict: Verdict,
}

#[derive(Serialize)]
pub struct Separated {
    pub stage: String,
    pub quantity: String,
    pub limit: f64,
    pub reference: f64,
    pub values: Vec<f64>,
    pub discrete_relative_errors: Vec<f64>,
    pub discrete_verdicts: Vec<Verdict>,
    /// Existing harness extrapolation. When every value lies on one side of
    /// the reference, this estimate equals the reference by algebra. It is
    /// recorded and is not the 5% decision.
    pub harness_extrapolation: Extrapolation,
    pub reference_free_order: Option<f64>,
    pub reference_free_estimate: Option<f64>,
    pub discretization_component: Option<f64>,
    pub continuum_relative_error: Option<f64>,
    pub continuum_verdict: Verdict,
    pub note: String,
}

#[derive(Serialize)]
pub struct ElementSeries {
    pub element: String,
    pub levels: Vec<Level>,
    pub crest: Vec<Separated>,
    pub field_relative_l2: Vec<Separated>,
    pub finest_mesh_five_percent: Verdict,
    pub continuum_five_percent_after_separation: Verdict,
}

#[derive(Serialize)]
pub struct Report {
    pub schema: String,
    pub profile: String,
    pub p6_complete: bool,
    pub griffiths_slope_five_percent_closed: bool,
    pub displacement_sample_relative_l2_limit: f64,
    pub maximum_displacement_error_over_reference_max: f64,
    pub depth_variation_limit: f64,
    pub original_2d_crest_m: [f64; 3],
    pub original_2d_theory_displacement_error_m: f64,
    pub documented_crest_matches_closed_form: bool,
    pub column_problem_only: bool,
    pub series: Vec<ElementSeries>,
    pub column_finest_mesh_five_percent: Verdict,
    pub column_continuum_five_percent_after_separation: Verdict,
    pub archived_slope: Vec<ArchivedSlope>,
    pub archived_slope_continuum_five_percent_verdict: Verdict,
    pub limitations: Vec<String>,
    pub user_command: String,
    pub slope_user_procedure: String,
}

#[derive(Debug, Deserialize)]
struct SeriesInput {
    schema: String,
    quantity: String,
    reference: f64,
    nested_h_refinement: bool,
    values: Vec<f64>,
}

#[derive(Serialize)]
pub struct SeriesJudgment {
    pub schema: String,
    pub quantity: String,
    pub nested_h_refinement: bool,
    pub p6_complete: bool,
    pub judgment: Separated,
}

pub fn judge_series(input: &str) -> Result<SeriesJudgment> {
    let series: SeriesInput = serde_json::from_str(input)?;
    if series.schema != SERIES_SCHEMA {
        return Err(Error::Input(
            "plane-strain series schema is not soil-fem-plane-strain-2d-series/1".into(),
        ));
    }
    if series.quantity.is_empty() || series.values.len() < 2 || !series.reference.is_finite() {
        return Err(Error::Input(
            "a plane-strain series needs a name, a finite reference and at least two values".into(),
        ));
    }
    if series.reference.abs() <= f64::MIN_POSITIVE {
        return Err(Error::Input(
            "the 5% ratio needs a non-zero original 2D displacement reference".into(),
        ));
    }
    if series.values.iter().any(|value| !value.is_finite()) {
        return Err(Error::Input(
            "plane-strain series values must be finite".into(),
        ));
    }
    let judgment = if series.nested_h_refinement {
        separate(&series.quantity, series.reference, &series.values, false)
    } else {
        separate_discrete_only(&series.quantity, series.reference, &series.values)
    };
    Ok(SeriesJudgment {
        schema: SCHEMA.into(),
        quantity: series.quantity,
        nested_h_refinement: series.nested_h_refinement,
        p6_complete: false,
        judgment,
    })
}

pub fn run(profile: Profile) -> Result<Report> {
    if !documented_crest_agrees_with_closed_form()? {
        return Err(Error::Input(
            "printed original 2D crest values disagree with the closed form beyond the archived 4.14e-13 m".into(),
        ));
    }
    let (tet4, tet10) = match profile {
        Profile::Ci => (vec![1, 2, 4], vec![1, 2, 4]),
        Profile::User => (vec![2, 4, 8], vec![1, 2, 4]),
    };
    let series = vec![
        element_series("tet4", false, &tet4)?,
        element_series("tet10", true, &tet10)?,
    ];
    let column_finest_mesh_five_percent =
        combine(series.iter().map(|item| item.finest_mesh_five_percent));
    let column_continuum_five_percent_after_separation = combine(
        series
            .iter()
            .map(|item| item.continuum_five_percent_after_separation),
    );
    let archived_slope: Vec<_> = ARCHIVED_SLOPE_MIDDLE_RELATIVE_L2
        .iter()
        .map(|(case, value)| ArchivedSlope {
            case: (*case).into(),
            middle_depth_displacement_sample_relative_l2: *value,
            five_percent_verdict: if *value <= DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT {
                Verdict::Pass
            } else {
                Verdict::Fail
            },
            in_plane_h_series: false,
            continuum_five_percent_verdict: Verdict::Indeterminate,
        })
        .collect();
    Ok(Report {
        schema: SCHEMA.into(),
        profile: match profile {
            Profile::Ci => "ci",
            Profile::User => "user",
        }
        .into(),
        p6_complete: false,
        griffiths_slope_five_percent_closed: false,
        displacement_sample_relative_l2_limit: DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT,
        maximum_displacement_error_over_reference_max:
            MAXIMUM_DISPLACEMENT_ERROR_OVER_REFERENCE_MAX,
        depth_variation_limit: DEPTH_DISPLACEMENT_VARIATION_LIMIT,
        original_2d_crest_m: ORIGINAL_2D_CREST_M,
        original_2d_theory_displacement_error_m: ORIGINAL_2D_THEORY_DISPLACEMENT_ERROR_M,
        documented_crest_matches_closed_form: true,
        column_problem_only: true,
        series,
        column_finest_mesh_five_percent,
        column_continuum_five_percent_after_separation,
        archived_slope,
        archived_slope_continuum_five_percent_verdict: Verdict::Indeterminate,
        limitations: limitations(),
        user_command: "cargo run --release -- plane-strain-2d user runs/plane-strain-2d-user.json"
            .into(),
        slope_user_procedure: slope_procedure(),
    })
}

fn limitations() -> Vec<String> {
    vec![
        "The 5% limit is the archived displacement_sample_relative_L2 of 0.05. It is not relaxed.".into(),
        "Printed original 2D crest values and the closed form they matched are not rewritten.".into(),
        "Archived Griffiths failures near 15.6% to 15.9% stay failures. Depth, load-increment and material-substep changes are not an in-plane h series, so Richardson is not applied to them and the slope continuum verdict stays indeterminate.".into(),
        "A column pass is this constrained MC column only. It does not close the slope, excavation meshes from 274 to 14909 elements, general 3D convergence, Fs, or P6.".into(),
        "The 5% decision uses a reference-free Richardson estimate: the order comes from the last two differences of the computed values, then that order is applied to the two finest values. spatial_convergence::extrapolate is recorded as well. When every sample lies on one side of the 2D value, that harness estimate equals the 2D value by algebra and is not the decision.".into(),
        "Depth-face variation uses the archived 1% companion and is reported separately. A Tet4 failure of that companion is not deleted and does not relax the limit.".into(),
    ]
}

fn slope_procedure() -> String {
    "The published 661-node 2D displacement file is not in this repository. Do not rerun Excel. On a workstation, expand griffiths_1999_coarse from pinned revision 440d246d9e32f5892c79539d4a7b797ad29442ee with interim/2026-10-30-Root635-full-audits-slope-reference/evidence-drivers/prepare-2d-slope-reference.py, build the plane-strain extrusions with prepare-3d-slope-match.py, and solve xy1-ny1 and xy2-ny2 at the same depth count. Do not mix ny changes into that pair. Pass one shared scalar displacement (the same material point or node on every in-plane mesh) and its original 2D value as a nested series: cargo run --release -- plane-strain-2d judge SERIES.json NEW_REPORT.json. Set nested_h_refinement true only for that in-plane pair. Leave it false for the archived depth and increment records.".into()
}

fn combine(verdicts: impl Iterator<Item = Verdict>) -> Verdict {
    let mut seen = false;
    let mut result = Verdict::Pass;
    for verdict in verdicts {
        seen = true;
        result = match (result, verdict) {
            (Verdict::Indeterminate, _) | (_, Verdict::Indeterminate) => Verdict::Indeterminate,
            (Verdict::Fail, _) | (_, Verdict::Fail) => Verdict::Fail,
            (Verdict::Pass, Verdict::Pass) => Verdict::Pass,
        };
    }
    if seen { result } else { Verdict::Indeterminate }
}

fn element_series(element: &str, tet10: bool, divisions: &[usize]) -> Result<ElementSeries> {
    let mut levels = Vec::new();
    for &division in divisions {
        levels.push(column_level(division, tet10)?);
    }
    let mut crest = Vec::new();
    let mut field_relative_l2 = Vec::new();
    for stage in 0..3 {
        let crests: Vec<_> = levels
            .iter()
            .map(|level| level.stages[stage].crest_mean_uz_m)
            .collect();
        let samples: Vec<_> = levels
            .iter()
            .map(|level| level.stages[stage].displacement_sample_relative_l2)
            .collect();
        crest.push(separate(
            &format!("{element} {} crest mean uz", STAGE_NAMES[stage]),
            ORIGINAL_2D_CREST_M[stage],
            &crests,
            false,
        ));
        field_relative_l2.push(separate(
            &format!(
                "{element} {} displacement sample relative L2",
                STAGE_NAMES[stage]
            ),
            0.,
            &samples,
            true,
        ));
    }
    let finest_mesh_five_percent = if levels.last().is_some_and(|level| {
        level
            .stages
            .iter()
            .all(|sample| sample.relative_l2_verdict == Verdict::Pass)
    }) {
        Verdict::Pass
    } else {
        Verdict::Fail
    };
    let continuum_five_percent_after_separation = combine(
        crest
            .iter()
            .chain(&field_relative_l2)
            .map(|item| item.continuum_verdict),
    );
    Ok(ElementSeries {
        element: element.into(),
        levels,
        crest,
        field_relative_l2,
        finest_mesh_five_percent,
        continuum_five_percent_after_separation,
    })
}

fn separate(quantity: &str, reference: f64, values: &[f64], already_relative: bool) -> Separated {
    let scale = if already_relative {
        1.
    } else {
        reference.abs()
    };
    let harness_reference = if already_relative { 0. } else { reference };
    let harness_extrapolation = spatial_convergence::extrapolate(
        quantity,
        harness_reference,
        values.to_vec(),
        1e-12,
        "Recorded with spatial_convergence::extrapolate. A one-sided approach makes this estimate equal the reference, so it is not the 5% decision.",
    );
    let discrete_relative_errors: Vec<_> = if already_relative {
        values.to_vec()
    } else {
        values
            .iter()
            .map(|value| (value - reference).abs() / scale)
            .collect()
    };
    let discrete_verdicts = discrete_relative_errors
        .iter()
        .copied()
        .map(|error| ratio_verdict(error, DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT))
        .collect();
    let at_roundoff = !discrete_relative_errors.is_empty()
        && discrete_relative_errors.iter().all(|error| *error <= 1e-8);
    let (reference_free_order, reference_free_estimate) = if at_roundoff {
        (None, None)
    } else {
        reference_free_richardson(values)
    };
    let (discretization_component, continuum_relative_error, continuum_verdict, note) = if let (
        Some(order),
        Some(estimate),
    ) =
        (reference_free_order, reference_free_estimate)
    {
        let continuum = if already_relative {
            estimate.abs()
        } else {
            (estimate - reference).abs() / scale
        };
        let component = values.last().map(|fine| (fine - estimate).abs() / scale);
        (
            component,
            Some(continuum),
            ratio_verdict(continuum, DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT),
            format!(
                "Reference-free order {order:.6} from the last two mesh differences. The 5% decision uses this estimate."
            ),
        )
    } else if !discrete_relative_errors.is_empty()
        && discrete_relative_errors
            .iter()
            .all(|error| *error <= DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT)
        && values
            .windows(2)
            .all(|pair| (pair[0] - pair[1]).abs() / scale <= 1e-8)
    {
        let fine = *values.last().expect("values");
        let continuum = if already_relative {
            fine.abs()
        } else {
            (fine - reference).abs() / scale
        };
        (
                Some(0.),
                Some(continuum),
                ratio_verdict(continuum, DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT),
                "Successive changes sit at the solver floor, so no order is identified. The estimate is the finest value.".into(),
            )
    } else {
        (
                None,
                None,
                Verdict::Indeterminate,
                "Mesh differences do not identify an order in [0.5, 8]. No continuum 5% verdict is issued.".into(),
            )
    };
    Separated {
        stage: quantity.into(),
        quantity: quantity.into(),
        limit: DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT,
        reference,
        values: values.to_vec(),
        discrete_relative_errors,
        discrete_verdicts,
        harness_extrapolation,
        reference_free_order,
        reference_free_estimate,
        discretization_component,
        continuum_relative_error,
        continuum_verdict,
        note,
    }
}

fn separate_discrete_only(quantity: &str, reference: f64, values: &[f64]) -> Separated {
    let scale = reference.abs();
    let errors: Vec<_> = values
        .iter()
        .map(|value| (value - reference).abs())
        .collect();
    let harness_extrapolation = Extrapolation {
        quantity: quantity.into(),
        reference,
        values: values.to_vec(),
        errors,
        observed_orders: vec![None; values.len().saturating_sub(1)],
        richardson_estimate: None,
        richardson_order: None,
        estimate_minus_reference: None,
        order_observable: false,
        note: "This series is not a nested in-plane h refinement. Richardson is not applied."
            .into(),
    };
    let discrete_relative_errors: Vec<_> = values
        .iter()
        .map(|value| (value - reference).abs() / scale)
        .collect();
    let discrete_verdicts = discrete_relative_errors
        .iter()
        .copied()
        .map(|error| ratio_verdict(error, DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT))
        .collect();
    Separated {
        stage: quantity.into(),
        quantity: quantity.into(),
        limit: DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT,
        reference,
        values: values.to_vec(),
        discrete_relative_errors,
        discrete_verdicts,
        harness_extrapolation,
        reference_free_order: None,
        reference_free_estimate: None,
        discretization_component: None,
        continuum_relative_error: None,
        continuum_verdict: Verdict::Indeterminate,
        note: "Not a nested in-plane h series. Continuum 5% stays indeterminate. Discrete ratios are still recorded.".into(),
    }
}

/// Order from the last two solution differences, then Richardson on the two
/// finest values. The original 2D number is not used to form the order.
fn reference_free_richardson(values: &[f64]) -> (Option<f64>, Option<f64>) {
    if values.len() < 3 {
        return (None, None);
    }
    let count = values.len();
    let coarse = values[count - 3];
    let mid = values[count - 2];
    let fine = values[count - 1];
    let first = (coarse - mid).abs();
    let second = (mid - fine).abs();
    if !(first > 0. && second > 0. && second < 0.7 * first) {
        return (None, None);
    }
    let order = (first / second).ln() / 2_f64.ln();
    if !(0.5..=8.).contains(&order) {
        return (None, None);
    }
    let estimate = spatial_convergence::richardson(mid, fine, 2., order);
    (estimate.map(|_| order), estimate)
}

fn ratio_verdict(error: f64, limit: f64) -> Verdict {
    if error.is_finite() && error <= limit {
        Verdict::Pass
    } else {
        Verdict::Fail
    }
}

fn column_level(divisions: usize, tet10: bool) -> Result<Level> {
    let mut mesh = column_mesh(divisions)?;
    if tet10 {
        mesh = mesh.elevate_straight()?;
        mesh.quadrature_order = Some(5);
        mesh.validate()?;
    }
    let (ground, case) = column_inputs()?;
    let solved = stages::solve(&mesh, ground, case)?;
    if !solved.complete {
        return Err(Error::Solve(format!(
            "plane-strain column divisions {divisions} did not finish: {}",
            solved.failure_reason.unwrap_or_else(|| "incomplete".into())
        )));
    }
    if solved.frames.len() != 3 {
        return Err(Error::Solve(
            "plane-strain column did not return gravity, load and unload frames".into(),
        ));
    }
    let expected = ["initial_equilibrated", "apply", "unload"];
    let mut stages = Vec::new();
    for (stage, frame) in solved.frames.iter().enumerate() {
        if frame.name != expected[stage] {
            return Err(Error::Solve(format!(
                "unexpected column frame {} at stage {stage}",
                frame.name
            )));
        }
        stages.push(sample_field(&mesh, &frame.state.displacement_m, stage)?);
    }
    let column = original_column();
    Ok(Level {
        divisions,
        element: if tet10 { "tet10" } else { "tet4" }.into(),
        h_m: column.width_m / divisions as f64,
        nodes: mesh.nodes.len(),
        elements: mesh.tetrahedra.len(),
        stages,
    })
}

fn sample_field(mesh: &Mesh, displacement: &[[f64; 3]], stage: usize) -> Result<FieldSample> {
    if displacement.len() != mesh.nodes.len() {
        return Err(Error::Solve(
            "plane-strain displacement length does not match the mesh".into(),
        ));
    }
    let mut numerator = 0_f64;
    let mut denominator = 0_f64;
    let mut max_error = 0_f64;
    let mut max_reference = 0_f64;
    let mut max_depth = 0_f64;
    let mut max_lateral = 0_f64;
    let mut crest_sum = 0_f64;
    let mut crest_count = 0_f64;
    let height = original_column().height_m;
    for (node, point) in mesh.nodes.iter().enumerate() {
        let reference = original_2d_vertical_displacement_m(point[2], stage)?;
        let computed = displacement[node];
        max_lateral = max_lateral.max(computed[0].abs());
        max_depth = max_depth.max(computed[1].abs());
        let difference = [computed[0], computed[2] - reference];
        let error = difference[0].hypot(difference[1]);
        let reference_norm = reference.abs();
        numerator += error * error;
        denominator += reference_norm * reference_norm;
        max_error = max_error.max(error);
        max_reference = max_reference.max(reference_norm);
        if (point[2] - height).abs() <= 1e-12 {
            crest_sum += computed[2];
            crest_count += 1.;
        }
    }
    if crest_count == 0. || denominator <= 0. || max_reference <= 0. {
        return Err(Error::Solve(
            "plane-strain column has no crest or no reference displacement".into(),
        ));
    }
    let relative_l2 = (numerator / denominator).sqrt();
    let maximum_ratio = max_error / max_reference;
    let depth_variation = depth_variation(mesh, displacement) / max_reference;
    Ok(FieldSample {
        displacement_sample_relative_l2: relative_l2,
        maximum_displacement_error_over_reference_max: maximum_ratio,
        depth_displacement_variation_over_reference_max: depth_variation,
        maximum_abs_depth_displacement_m: max_depth,
        maximum_abs_lateral_displacement_m: max_lateral,
        crest_mean_uz_m: crest_sum / crest_count,
        crest_reference_m: ORIGINAL_2D_CREST_M[stage],
        relative_l2_verdict: ratio_verdict(relative_l2, DISPLACEMENT_SAMPLE_RELATIVE_L2_LIMIT),
        maximum_ratio_verdict: ratio_verdict(
            maximum_ratio,
            MAXIMUM_DISPLACEMENT_ERROR_OVER_REFERENCE_MAX,
        ),
        depth_variation_verdict: ratio_verdict(depth_variation, DEPTH_DISPLACEMENT_VARIATION_LIMIT),
    })
}

fn depth_variation(mesh: &Mesh, displacement: &[[f64; 3]]) -> f64 {
    let mut groups: BTreeMap<(i64, i64), Vec<usize>> = BTreeMap::new();
    for (index, point) in mesh.nodes.iter().enumerate() {
        let key = (
            (point[0] * 1e9).round() as i64,
            (point[2] * 1e9).round() as i64,
        );
        groups.entry(key).or_default().push(index);
    }
    let mut worst = 0_f64;
    for indices in groups.values() {
        for left in indices {
            for right in indices {
                let a = displacement[*left];
                let b = displacement[*right];
                worst = worst.max((a[0] - b[0]).abs()).max((a[2] - b[2]).abs());
            }
        }
    }
    worst
}

fn column_inputs() -> Result<(Ground, StageCase)> {
    let column = original_column();
    let material = mc::Parameters {
        young_kpa: column.young_kpa,
        poisson: column.poisson,
        cohesion_kpa: column.cohesion_kpa,
        friction_deg: column.friction_deg,
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
    states.insert("Column".into(), zero);
    let ground = Ground {
        schema: "soil-fem-ground/1".into(),
        units: "m-kN-kPa".into(),
        initialization: Initialization {
            method: "equilibrate".into(),
            surface_surcharge_kpa: 0.,
            water_level_m: None,
            water_unit_weight_kn_m3: 9.80665,
            relative_equilibrium_tolerance: 1e-8,
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
            name: "original 2D MC column".into(),
            unit_weight_kn_m3: column.unit_weight_kn_m3,
            model: soil::Model::MohrCoulomb(material),
        }],
        construction_material_ids: vec![],
        layers: vec![],
        bodies: vec![geostatic::Body {
            body: "Column".into(),
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
            initial_increment: 1. / 16.,
            minimum_increment: 1. / 1024.,
            material_maximum_substep: 1e-4,
            mcc_integration: None,
            tangent_perturbation: 1e-7,
            maximum_free_dofs: 100_000,
            linear_solver: LinearSolver::SparseGmresIlu0 {
                options: GmresOptions::default(),
            },
        },
        boundary: Some(column_boundary()),
        monitors: BTreeMap::new(),
        monitor_sets: BTreeMap::new(),
        pressure_on_original_surface_only: true,
        stages: vec![
            Stage::Pressure {
                name: "apply".into(),
                target_pressure_kpa: column.pressure_kpa,
                initial_increment: Some(1. / 16.),
            },
            Stage::Pressure {
                name: "unload".into(),
                target_pressure_kpa: 0.,
                initial_increment: Some(1. / 16.),
            },
        ],
        initially_inactive_bodies: vec![],
        verification_reject_once_stages: vec![],
    };
    Ok((ground, case))
}

fn column_boundary() -> Conditions {
    let fixed = |set: &str, axis: usize| {
        let mut direction = [0.; 3];
        direction[axis] = 1.;
        Constraint {
            node_set: set.into(),
            direction,
            value_m: 0.,
            gradient: [0.; 3],
            reference_point_m: [0.; 3],
        }
    };
    Conditions {
        constraints: vec![fixed("all", 0), fixed("all", 1), fixed("bottom", 2)],
        nodal_loads: vec![],
        surface_loads: vec![],
    }
}

fn column_mesh(divisions: usize) -> Result<Mesh> {
    if !(1..=32).contains(&divisions) {
        return Err(Error::Input(
            "plane-strain column divisions must be 1..=32".into(),
        ));
    }
    let column = original_column();
    let (nx, ny, nz) = (divisions, 1, divisions);
    let size = [column.width_m, column.depth_m, column.height_m];
    let mut nodes = Vec::new();
    for iz in 0..=nz {
        for iy in 0..=ny {
            for ix in 0..=nx {
                nodes.push([
                    size[0] * ix as f64 / nx as f64,
                    size[1] * iy as f64 / ny as f64,
                    size[2] * iz as f64 / nz as f64,
                ]);
            }
        }
    }
    let id = |ix: usize, iy: usize, iz: usize| ix + (nx + 1) * (iy + (ny + 1) * iz);
    let mut tetrahedra = Vec::new();
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
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
    let count = tetrahedra.len();
    let mut node_sets = BTreeMap::new();
    node_sets.insert("all".into(), (0..nodes.len()).collect());
    node_sets.insert(
        "bottom".into(),
        nodes
            .iter()
            .enumerate()
            .filter_map(|(index, point)| (point[2].abs() <= 1e-12).then_some(index))
            .collect(),
    );
    let hash = format!("{:x}", Sha256::digest(b"plane-strain-original-2d-column"));
    let volume = size[0] * size[1] * size[2];
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
        material_ids: vec![1; count],
        node_sets,
        face_sets: BTreeMap::new(),
        regions: vec![Region {
            id: 1,
            body: "Column".into(),
            material_id: 1,
            bounds_m: [[0.; 3], size],
            volume_m3: volume,
            source_sha256: hash,
        }],
        region_ids: vec![1; count],
        excavation_sets: BTreeMap::new(),
    };
    let measured = mesh.validate()?;
    if (measured - volume).abs() > 1e-9 {
        return Err(Error::Input(format!(
            "plane-strain column volume {measured} is not {volume}"
        )));
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn documented_crest_stays_inside_the_archived_theory_band() -> Result<()> {
        assert!(documented_crest_agrees_with_closed_form()?);
        let height = original_column().height_m;
        let gravity = original_2d_vertical_displacement_m(height, 0)?;
        assert!((gravity + original_column().height_m / 600.).abs() < 1e-15);
        Ok(())
    }

    #[test]
    fn nested_monitor_separates_a_coarse_miss_from_the_continuum() {
        let reference = ORIGINAL_2D_CREST_M[1];
        let values = [reference + 4e-4, reference + 1e-4, reference + 2.5e-5];
        let judgment = separate("synthetic loaded crest", reference, &values, false);
        assert_eq!(judgment.discrete_verdicts[0], Verdict::Fail);
        assert_eq!(judgment.discrete_verdicts[2], Verdict::Pass);
        assert_eq!(judgment.continuum_verdict, Verdict::Pass);
        assert!(judgment.continuum_relative_error.unwrap() < 1e-9);
        assert!(judgment.discretization_component.unwrap() > 0.);
        assert_eq!(judgment.limit.to_bits(), 0.05_f64.to_bits());
    }

    #[test]
    fn non_nested_series_stays_indeterminate_and_keeps_discrete_failures() {
        let input = r#"{
            "schema": "soil-fem-plane-strain-2d-series/1",
            "quantity": "archived slope middle depth is not an h series",
            "reference": -1.0,
            "nested_h_refinement": false,
            "values": [-0.8441800129308275, -0.8411430807048378]
        }"#;
        let judged = judge_series(input).expect("series");
        assert!(!judged.p6_complete);
        assert!(!judged.nested_h_refinement);
        assert_eq!(judged.judgment.continuum_verdict, Verdict::Indeterminate);
        assert!(
            judged
                .judgment
                .harness_extrapolation
                .richardson_estimate
                .is_none()
        );
        assert!(
            judged
                .judgment
                .discrete_verdicts
                .iter()
                .all(|verdict| *verdict == Verdict::Fail)
        );
    }

    #[test]
    fn plateau_away_from_the_2d_value_is_not_called_a_pass() {
        let reference = ORIGINAL_2D_CREST_M[0];
        let values = [reference + 8e-5, reference + 7.2e-5, reference + 7.1e-5];
        let judgment = separate("plateau", reference, &values, false);
        assert_eq!(judgment.continuum_verdict, Verdict::Fail);
        assert!(judgment.continuum_relative_error.unwrap() > 0.05);
        assert!(
            judgment
                .discrete_verdicts
                .iter()
                .all(|verdict| *verdict == Verdict::Fail)
        );
    }

    #[test]
    fn ci_column_records_five_percent_without_closing_the_slope_or_p6() -> Result<()> {
        let report = run(Profile::Ci)?;
        let text = serde_json::to_string_pretty(&report)?;
        assert!(!report.p6_complete, "{text}");
        assert!(!report.griffiths_slope_five_percent_closed);
        assert_eq!(
            report.displacement_sample_relative_l2_limit.to_bits(),
            0.05_f64.to_bits()
        );
        assert!(report.archived_slope.iter().all(|row| {
            row.five_percent_verdict == Verdict::Fail
                && row.continuum_five_percent_verdict == Verdict::Indeterminate
                && !row.in_plane_h_series
                && row.middle_depth_displacement_sample_relative_l2 > 0.05
        }));
        assert_eq!(report.series.len(), 2);
        let tet4 = report
            .series
            .iter()
            .find(|series| series.element == "tet4")
            .unwrap();
        let tet10 = report
            .series
            .iter()
            .find(|series| series.element == "tet10")
            .unwrap();
        assert_eq!(tet4.finest_mesh_five_percent, Verdict::Fail);
        assert_eq!(tet4.continuum_five_percent_after_separation, Verdict::Pass);
        assert_eq!(tet10.finest_mesh_five_percent, Verdict::Pass);
        assert_eq!(tet10.continuum_five_percent_after_separation, Verdict::Pass);
        assert_eq!(report.column_finest_mesh_five_percent, Verdict::Fail);
        assert_eq!(
            report.column_continuum_five_percent_after_separation,
            Verdict::Pass
        );
        let gravity = &tet4
            .levels
            .iter()
            .find(|level| level.divisions == 4)
            .unwrap()
            .stages[0];
        assert!((gravity.displacement_sample_relative_l2 - 0.06828189917309707).abs() < 1e-12);
        assert_eq!(gravity.relative_l2_verdict, Verdict::Fail);
        assert_eq!(gravity.depth_variation_verdict, Verdict::Fail);
        let separated = tet4
            .field_relative_l2
            .iter()
            .find(|item| item.stage.contains("gravity"))
            .unwrap();
        let continuum = separated.continuum_relative_error.unwrap();
        assert!((continuum - 0.010335840557563397).abs() < 1e-12, "{text}");
        assert!(continuum < 0.05);
        for level in tet4.levels.iter().chain(&tet10.levels) {
            for sample in &level.stages {
                assert!(sample.maximum_abs_depth_displacement_m <= 1e-12);
                assert!(sample.maximum_abs_lateral_displacement_m <= 1e-12);
            }
        }
        let tet10_fine = tet10.levels.last().unwrap();
        assert!(tet10_fine.divisions == 4);
        assert!(tet10_fine.stages.iter().all(|sample| {
            sample.relative_l2_verdict == Verdict::Pass
                && sample.displacement_sample_relative_l2 < 1e-4
        }));
        Ok(())
    }

    #[test]
    fn archived_slope_failure_files_are_unchanged() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let files = [
            (
                "xy1-ny1",
                "interim/2026-10-30-Root635-full-audits-slope-reference/slope/coarse-field-FAILED.json",
            ),
            (
                "xy1-ny1",
                "interim/2026-10-30-Root635-CI-and-refinement-controls/slope/xy1-ny1/field-acceptance-FAILED.json",
            ),
            (
                "xy1-ny2",
                "interim/2026-10-30-Root635-CI-and-refinement-controls/slope/xy1-ny2/field-acceptance-FAILED.json",
            ),
            (
                "xy1-ny1-global64",
                "interim/2026-10-30-Root635-CI-and-refinement-controls/slope/xy1-ny1-global64/field-acceptance-FAILED.json",
            ),
            (
                "xy1-ny1-local-half",
                "interim/2026-10-30-async-UI-and-independent-SRM-bound/slope/material-half1200/field-acceptance-FAILED.json",
            ),
        ];
        for (case, relative) in files {
            let value: serde_json::Value = serde_json::from_str(
                &std::fs::read_to_string(root.join(relative))
                    .unwrap_or_else(|error| panic!("{relative}: {error}")),
            )
            .unwrap();
            assert_eq!(value["passed"], false, "{relative}");
            assert_eq!(
                value["gates"]["all_661_displacements_each_depth"], false,
                "{relative}"
            );
            assert_eq!(
                value["original_predeclared_limits"]["displacement_sample_relative_L2"],
                0.05
            );
            let middle = value["metrics"]["displacement_sample_relative_L2_per_depth"]["0.5"]
                .as_f64()
                .unwrap();
            let archived = ARCHIVED_SLOPE_MIDDLE_RELATIVE_L2
                .iter()
                .find(|(name, _)| *name == case)
                .unwrap()
                .1;
            assert!(
                (middle - archived).abs() <= 1e-15,
                "{case} {middle} {archived}"
            );
            assert!(middle > 0.05);
        }
    }
}
