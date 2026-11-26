        if coarse.state.stage_index == 3 && !audit.reject
            && let Some(path) = std::env::var_os("P9_MATCHED") {
            type Targets = Vec<(f64, f64, Vec<f64>)>;
            static TABLE: std::sync::OnceLock<Targets> = std::sync::OnceLock::new();
            let rows = TABLE.get_or_init(|| {
                let f = std::fs::File::open(path).unwrap();
                let plan: serde_json::Value = serde_json::from_reader(f).unwrap();
                let level = std::env::var("P9_LEVEL").unwrap();
                serde_json::from_value(plan["levels"][level].clone()).unwrap()
            });
            if let Some((_, _, times)) = rows.iter().find(|(l,r,_)| *l==start && *r==load.progress) {
                for &time in times {
                    let mid = time <= midpoint_progress;
                    let solved;
                    let (j,residual) = if time == midpoint_progress {
                        (b,half.residual)
                    } else if time == load.progress {
                        (d,fine.residual)
                    } else {
                        let f=(time-start)/h;
                        let mut target=load.clone();
                        target.progress=time;
                        target.boundary_factor=old.boundary_factor+(load.boundary_factor-old.boundary_factor)*f;
                        target.pressure=old.surface_pressure_kpa+(load.pressure-old.surface_pressure_kpa)*f;
                        for ((v,a),b) in target.load.iter_mut().zip(&old.applied_nodal_load_kn).zip(&load.load) {
                            *v=a+(b-a)*f;
                        }
                        for ((v,a),b) in target.release.iter_mut().zip(old_release).zip(&load.release) {
                            *v=a+(b-a)*f;
                        }
                        solved=endpoint(if mid { old } else { &half.state },target,true)?.0;
                        (solved.state.joints.as_ref().unwrap(),solved.residual)
                    };
                    let points=[3,5,13].map(|i| {
                        let p=&j.points[i];
                        serde_json::json!([p.status,p.plastic_slip_m])
                    });
                    eprintln!("{}",serde_json::json!({"matched":time,"parent":[start,midpoint_progress,load.progress],
                                              "first_half":mid,"force":residual,"points":points}));
                }
            }
        }
