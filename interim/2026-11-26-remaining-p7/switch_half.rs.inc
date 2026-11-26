// Only the P7 experiment sets this opt-in. Preserve the original full/two-half audit.
if p7 && coarse.state.stage_index == 3 && !audit.reject {
    let switches = [a.points.iter().zip(&b.points).any(|(u,v)| u.status != v.status),
                    b.points.iter().zip(&d.points).any(|(u,v)| u.status != v.status)];
    let mut extra_reject = false;
    for (side, switched) in switches.into_iter().enumerate() {
        if !switched { continue; }
        // Both halves have equal length. Reserve half the fixed slip and work budgets
        // for each (including when both halves contain a transition).
        let f = (side as f64 + 0.5) / 2.;
        let mut q = load.clone();
        q.progress = start + h*f;
        q.boundary_factor = old.boundary_factor + (load.boundary_factor-old.boundary_factor)*f;
        q.pressure = old.surface_pressure_kpa + (load.pressure-old.surface_pressure_kpa)*f;
        for ((v,a),b) in q.load.iter_mut().zip(&old.applied_nodal_load_kn).zip(&load.load) {
            *v = a+(b-a)*f;
        }
        for ((v,a),b) in q.release.iter_mut().zip(old_release).zip(&load.release) {
            *v = a+(b-a)*f;
        }
        let (parent, end, reference) = if side==0 {
            (old, mid_load.clone(), b)
        } else {
            (&half.state, load.clone(), d)
        };
        let quarter = endpoint(parent,q,true)?.0;
        let refined = endpoint(&quarter.state,end,true)?.0;
        let mut slip: f64 = 0.;
        let mut work: f64 = 0.;
        for (u,v) in refined.state.joints.as_ref().unwrap().points.iter().zip(&reference.points) {
            for i in 0..3 {
                slip = slip.max((u.plastic_slip_m[i]-v.plastic_slip_m[i]).abs());
            }
            work = work.max((u.plastic_dissipation_kpa_m-v.plastic_dissipation_kpa_m).abs());
        }
        let outside = slip > slip_budget_m*0.5 || work > dissipation_budget_kpa_m*0.5;
        eprintln!("{}",serde_json::json!({"half":side,"from":start,"to":start+h,
            "slip":slip,"work":work,"slip_budget":slip_budget_m*0.5,
            "work_budget":dissipation_budget_kpa_m*0.5,"outside":outside}));
        extra_reject |= outside;
    }
    if extra_reject {
        audit.reject = true;
        audit.proposed_increment = Some(h*0.5);
    }
}
