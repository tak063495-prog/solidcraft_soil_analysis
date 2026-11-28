// Diagnostic only. The production audit and returned Trial are already determined.
let windows: Vec<[f64; 2]> = serde_json::from_str(
    &std::env::var("P6_WINDOWS").expect("diagnostic windows")
).unwrap();
if windows.contains(&[start, start + h]) {
    let saved = serde_json::to_vec(old).unwrap();
    let restored: Committed = serde_json::from_slice(&saved).unwrap();
    assert_eq!(saved, serde_json::to_vec(&restored).unwrap());
    let at = |f: f64| {
        let mut x = load.clone();
        if f == 1. { return x; }
        x.progress = start + h * f;
        x.boundary_factor = old.boundary_factor + (load.boundary_factor-old.boundary_factor)*f;
        x.pressure = old.surface_pressure_kpa + (load.pressure-old.surface_pressure_kpa)*f;
        for ((v,a),b) in x.load.iter_mut().zip(&old.applied_nodal_load_kn).zip(&load.load) {
            *v = a + (b-a)*f;
        }
        for ((v,a),b) in x.release.iter_mut().zip(old_release).zip(&load.release) {
            *v = a + (b-a)*f;
        }
        x
    };
    let mut calls = 0;
    let mut solve = |p: &Committed, f: f64| -> Result<Trial> {
        let (t, corr, _) = endpoint(p, at(f))?;
        calls += 1 + usize::from(corr.is_some());
        Ok(t)
    };
    // Early allocation: 1/4, 1/2, 1. Late allocation: 1/2, 3/4, 1.
    // Both have three sequential equilibria. Recompute both from the serialized parent.
    let q = solve(&restored, 0.25)?;
    let mid = solve(&q.state, 0.5)?;
    let early = solve(&mid.state, 1.)?;
    let late_half = solve(&restored, 0.5)?;
    assert_eq!(serde_json::to_vec(&late_half.state).unwrap(), serde_json::to_vec(&half.state).unwrap());
    let q3 = solve(&late_half.state, 0.75)?;
    let late = solve(&q3.state, 1.)?;
    let snapshot = |s: &Committed| {
        let pts = &s.joints.as_ref().unwrap().points;
        serde_json::json!({
            "slip":pts.iter().map(|p|p.plastic_slip_m).collect::<Vec<_>>(),
            "work":pts.iter().map(|p|p.plastic_dissipation_kpa_m).collect::<Vec<_>>(),
            "status":pts.iter().map(|p|&p.status).collect::<Vec<_>>()
        })
    };
    eprintln!("{}", serde_json::json!({"allocation":true, "start":start,"end":start+h,
        "saved_parent":restored,"parent":snapshot(old),"coarse":snapshot(&coarse.state),"half":snapshot(&half.state),
        "fine":snapshot(&fine.state), "q1":snapshot(&q.state),"early_mid":snapshot(&mid.state),
        "early":snapshot(&early.state),"q3":snapshot(&q3.state),"late":snapshot(&late.state),
        "extra_solves":calls,"event":audit.friction_event_target, "accepted":!audit.reject,
        "slip_budget":slip_budget_m,"work_budget":dissipation_budget_kpa_m}));
}
