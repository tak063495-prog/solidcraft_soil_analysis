//! Tet4/Tet10 soil mechanics, geostatic initialization and staged MCC/MC FEM.
pub mod anchors;
pub mod beams;
pub mod benchmark;
pub mod boundary;
pub mod checkpoint;
pub mod darcy;
pub mod elastic;
pub mod elements;
pub mod embedding;
pub mod epoch_srm_control;
pub mod epoch_srm_gate;
pub mod epoch_srm_gate_job;
pub mod epoch_srm_job;
pub mod epoch_srm_search;
pub mod geostatic;
pub mod ground;
pub mod hardening;
pub mod hydraulic_transition;
pub mod hydrostatic;
pub mod joint_law;
pub mod joint_water;
pub mod joints;
pub mod linalg;
pub mod mc;
pub mod mc_epoch_job;
pub mod mc_epoch_path;
pub mod mcc;
pub mod mcc_split;
pub mod mechanics;
pub mod mesh;
pub mod mixed_flow;
pub mod mixed_mcc;
pub mod mixed_mechanics;
pub mod one_way;
pub mod one_way_construction;
pub mod plane_strain_2d;
pub mod plates;
pub mod replay;
pub mod retention;
pub mod richards;
pub mod richards_checkpoint;
pub mod rigidity;
pub mod soil;
pub mod solver;
pub mod spatial_convergence;
pub mod srm_bound;
pub mod srm_search;
pub mod srm_search_control;
pub mod srm_search_job;
pub mod srm_virtual_work;
pub mod stage_control;
pub mod stage_job;
pub mod stages;
pub mod surface_projection;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid input: {0}")]
    Input(String),
    #[error("numerical solve failed: {0}")]
    Solve(String),
    #[error("material update failed: {0}")]
    Material(String),
    #[error("runtime control failed: {0}")]
    Control(String),
    #[error(transparent)]
    Linear(#[from] linalg::LinearError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}
pub type Result<T> = std::result::Result<T, Error>;
