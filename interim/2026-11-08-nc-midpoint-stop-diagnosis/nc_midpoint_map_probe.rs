//! Read-only tangent audit of the existing midpoint increment map.
use serde::Deserialize;
use soil_fem_3d::{Result, Error, mcc, elastic::Tensor};
use std::{fs,path::Path};
#[derive(Deserialize)]
struct Request { parameters:mcc::Parameters, state:mcc::State, strain:[f64;6], partitions:usize }
fn rows(d:&Tensor)->Vec<Vec<f64>> {(0..6).map(|i|(0..6).map(|j|d[(i,j)]).collect()).collect()}
fn measure(r:&Request)->Result<serde_json::Value>{
    let primal=r.parameters.integrate_partition(&r.state,r.strain,r.partitions)?;
    let (implicit,d)=r.parameters.implicit_stress_tangent_partition(&r.state,r.strain,r.partitions)?;
    let mut fds=Vec::new();
    for h in [1e-8,1e-9] {
        let (_,f)=r.parameters.stress_tangent_partition(&r.state,r.strain,r.partitions,h)?;
        let mut branches=Vec::new();
        for j in 0..6 {let mut a=r.strain;let mut b=r.strain;a[j]+=h;b[j]-=h;
            branches.push([r.parameters.integrate_partition(&r.state,a,r.partitions)?.plastic_substeps,
                r.parameters.integrate_partition(&r.state,b,r.partitions)?.plastic_substeps]);}
        fds.push(serde_json::json!({"width":h,"tangent":rows(&f),
            "relative_to_implicit":(f-d).norm()/d.norm(),"plus_minus_plastic_substeps":branches}));
    }
    Ok(serde_json::json!({"primal":primal.state,"plastic_substeps":primal.plastic_substeps,
        "implicit_primal_exact":serde_json::to_value(&primal.state)?==serde_json::to_value(&implicit.state)?,
        "implicit_tangent":rows(&d),"finite_differences":fds}))
}
fn main()->Result<()>{
 let a:Vec<_>=std::env::args().collect();
 if a.len()!=3 || Path::new(&a[2]).exists(){return Err(Error::Input("nc_midpoint_map_probe INPUT NEW_OUTPUT".into()));}
 let rs:Vec<Request>=serde_json::from_slice(&fs::read(&a[1])?)?;
 let out:Vec<_>=rs.iter().map(|r|match measure(r){Ok(v)=>v,Err(e)=>serde_json::json!({"error":e.to_string()})}).collect();
 fs::write(&a[2],serde_json::to_vec(&out)?)?;Ok(())
}
