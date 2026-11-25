            if self.mesh.tetrahedra.len() == 12 && stage == 5 && old.stage_index == stage
                && old.stage_progress >= 926. / 1024. && progress <= 927. / 1024.
                && relative <= settings.relative_force_tolerance
            {
                use std::io::Write;
                let mut f = std::fs::OpenOptions::new().create(true).append(true).open("c01-paths.jsonl")?;
                let row = serde_json::json!({
                    "from":old.stage_progress, "to":progress, "factor":boundary_factor,
                    "parent_u":old.displacement_m, "parent":old.joints.as_ref().unwrap().points,
                    "u":u, "end":eval.joints.as_ref().unwrap().state.points,
                    "counts":partitions[self.point_count..], "residual":relative
                });
                serde_json::to_writer(&mut f, &row)?;
                writeln!(f)?;
            }
            if self.mesh.tetrahedra.len() == 96 && stage > 0 && iteration > 0
                && [NC_PAIRS].contains(&(old.stage_progress,progress))
                && (relative <= settings.relative_force_tolerance || old.stage_progress == 0.5647809505462646)
            {
                let mut points = Vec::new();
                for (ei,e) in self.elements.iter().enumerate().filter(|(i,_)| active[*i]) {
                    let mat = &self.material_models[self.material_index(&old.element_material_ids,ei)?];
                    let du = DVector::from_fn(e.dofs.len(), |i,_| {
                        let g = e.dofs[i];
                        u[g/3][g%3] - old.displacement_m[g/3][g%3]
                    });
                    for (qi,q) in e.points.iter().enumerate() {
                        let id = e.point_start + qi;
                        let de = -&q.b * &du;
                        let strain = std::array::from_fn(|i| de[i]);
                        let state = &old.material_states[id];
                        let steps = partitions[id];
                        let (forward,c) = mat.stress_tangent_partition(state,strain,steps,settings.tangent_perturbation)?;
                        // Diagnostic Jacobian on the branch selected by the forward increment.
                        // If just one perturbation crosses yield, use the same-branch side.
                        let h = (de.amax()*0.01).clamp(1e-13,settings.tangent_perturbation);
                        let mut narrow = crate::elastic::Tensor::zeros();
                        let mut single = 0;
                        let mut unresolved = 0;
                        for j in 0..6 {
                            let mut plus = strain;
                            let mut minus = strain;
                            plus[j] += h;
                            minus[j] -= h;
                            let p = mat.integrate_partition(state,plus,steps)?;
                            let m = mat.integrate_partition(state,minus,steps)?;
                            let (sp,sm,width) = match (p.plastic_substeps==forward.plastic_substeps,
                                                       m.plastic_substeps==forward.plastic_substeps) {
                                (true,true) => (p.state.stress(),m.state.stress(),2.*h),
                                (true,false) => { single+=1; (p.state.stress(),forward.state.stress(),h) }
                                (false,true) => { single+=1; (forward.state.stress(),m.state.stress(),h) }
                                (false,false) => { unresolved+=1; (p.state.stress(),m.state.stress(),2.*h) }
                            };
                            for i in 0..6 { narrow[(i,j)] = (sp[i]-sm[i])/width; }
                        }
                        let reverse = mat.integrate_partition(state,strain.map(|x|-x),steps)?;
                        points.push(serde_json::json!({
                            "id":id,"steps":steps,"positive":forward.plastic_substeps,"negative":reverse.plastic_substeps,
                            "width":h,"single":single,"unresolved":unresolved,"strain":strain,"c":c.as_slice(),"narrow":narrow.as_slice()
                        }));
                    }
                }
                let record = serde_json::json!({
                    "from":old.stage_progress,"to":progress,"iteration":iteration,"relative":relative,
                    "parent_hash":crate::checkpoint::hash_json(old)?, "target_load":target,
                    "old_load":old.applied_nodal_load_kn, "u":u,"old_u":old.displacement_m,
                    "points":points, "rows":dofs.rows
                });
                let name = format!("nc-{:016x}-{:016x}-{iteration}.json",old.stage_progress.to_bits(),progress.to_bits());
                serde_json::to_writer(std::io::BufWriter::new(std::fs::File::create(name)?),&record)?;
            }
