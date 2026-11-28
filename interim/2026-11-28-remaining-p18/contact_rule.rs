            let next = if p17 && (0.62..0.6875).contains(&progress)
                && std::env::var_os("P18_OBSERVE").is_some() {
                use crate::joint_law::Status;
                use nalgebra::Vector3;
                let pts = &committed.joints.as_ref().unwrap().points;
                let jumps: Vec<_> = pts.iter().map(|p| p.elastic_jump()).collect();
                if p18_latest.as_ref().is_none_or(|(t,_)| *t != progress) {
                    p18_prior = p18_latest.take();
                    p18_latest = Some((progress,jumps.clone()));
                }
                let nominal = next;
                let mut proposal = nominal;
                if let Some((prev,old)) = &p18_prior && progress >= 0.625 {
                    let joints = &case.joint_system.as_ref().unwrap().joints;
                    assert_eq!(joints.len(),1);
                    let mat = &joints[0].material;
                    let h = nominal - progress;
                    let base = stage.step(&case.settings);
                    let mut forecasts = Vec::new();
                    for (i,p) in pts.iter().enumerate() {
                        let n = Vector3::from(p.normal);
                        let dj = (jumps[i]-old[i])*(h/(progress-prev));
                        let gap = n.dot(&jumps[i]);
                        let dg = n.dot(&dj);
                        let tau = Vector3::from(p.traction_kpa);
                        let shear = tau - n*tau.dot(&n);
                        let ds = mat.tangential_stiffness_kpa_m*(dj-n*dg);
                        let cap = mat.friction_coefficient*mat.normal_stiffness_kpa_m;
                        let mut add = |frac: f64, kind: &str, guard: f64| {
                            let offset = h*frac;
                            if offset > guard && h-offset > guard {
                                forecasts.push((i,kind.to_string(),progress+offset));
                                proposal = proposal.min(progress+offset);
                            }
                        };
                        if (gap > 0.) != (gap+dg > 0.) {
                            let guard = if p.status==Status::Slip && dg>0. {
                                case.settings.minimum_increment
                            } else {base/64.};
                            add(-gap/dg,"gap",guard);
                        }
                        if p.status==Status::Stick && gap<=0. {
                            let end = if gap+dg>0. {-gap/dg} else {1.};
                            let margin = |f:f64| (shear+ds*f).norm()+cap*(gap+f*dg);
                            if margin(0.)<0. && margin(end)>0. {
                                let (mut lo,mut hi)=(0.,end);
                                for _ in 0..48 {
                                    let m=(lo+hi)*0.5;
                                    if margin(m)>0. {hi=m;} else {lo=m;}
                                }
                                add((lo+hi)*0.5,"friction",base/64.);
                            }
                        }
                    }
                    let sample:Vec<_> = pts.iter().map(|p|
                        (&p.status,p.elastic_jump().as_slice().to_vec(),p.normal,p.traction_kpa)).collect();
                    eprintln!("{}",serde_json::json!({"p18_parent":progress,"attempt":count,
                        "previous":prev,"previous_jumps":old.iter().map(|v| v.as_slice()).collect::<Vec<_>>(),"points":sample,"nominal":nominal,
                        "forecasts":forecasts,"proposal":proposal}));
                }
                if std::env::var_os("P18_CONTROL").is_some() {proposal} else {nominal}
            } else {next};
