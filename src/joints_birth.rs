//! Explicit interface generations and exact recovery of pre-birth forces.
use super::*;

impl Model {
    pub fn validate_lifecycle(&self, state: &State) -> Result<()> {
        let enabled = self
            .configuration
            .joints
            .iter()
            .any(|d| d.birth_policy.is_some());
        if !enabled && state.lifecycle.is_empty() {
            return Ok(());
        }
        let split = state.lifecycle.first().is_some_and(|l| l.region.is_some());
        let count = if split {
            self.faces.len()
        } else {
            self.configuration.joints.len()
        };
        if state.lifecycle.len() != count || state.points.len() != self.points.len() {
            return Err(Error::Input("JOINT BIRTH lifecycle layout mismatch".into()));
        }
        for (slot, l) in state.lifecycle.iter().enumerate() {
            let di = if split {
                self.points[self.faces[slot].start].material
            } else {
                slot
            };
            let d = &self.configuration.joints[di];
            if l.region != split.then(|| self.face_region(slot)) {
                return Err(Error::Input("JOINT BIRTH face area/owners mismatch".into()));
            }
            if l.joint_id != d.id
                || l.policy != d.birth_policy
                || (l.active && l.tied)
                || (d.birth_policy.is_none()
                    && (!l.active || l.generation != 0 || !l.archived.is_empty()))
                || (d.birth_policy.is_some() && l.active && l.generation == 0)
            {
                return Err(Error::Input(
                    "JOINT BIRTH policy/generation/activity mismatch".into(),
                ));
            }
            let indices: Vec<_> = self
                .points
                .iter()
                .enumerate()
                .filter(|(_, p)| p.slot(state) == slot)
                .map(|(i, _)| i)
                .collect();
            let mut previous = 0;
            for a in &l.archived {
                if a.generation == 0
                    || a.generation != previous + 1
                    || a.generation > l.generation
                    || a.points.iter().map(|(i, _)| *i).collect::<Vec<_>>() != indices
                {
                    return Err(Error::Input(
                        "JOINT BIRTH archived generation mismatch".into(),
                    ));
                }
                if a.endpoint_work_kn_m != split.then(|| self.frozen_work(&a.points)) {
                    return Err(Error::Input("JOINT BIRTH frozen face work mismatch".into()));
                }
                for (i, s) in &a.points {
                    if s.normal != state.points[*i].normal || s.birth.is_none() {
                        return Err(Error::Input(
                            "JOINT BIRTH frozen point reference mismatch".into(),
                        ));
                    }
                    d.material.validate_state(s)?;
                }
                previous = a.generation;
            }
            if l.generation != previous + u64::from(l.active && d.birth_policy.is_some())
                || (!l.active
                    && l.archived
                        .last()
                        .is_some_and(|a| a.points.iter().any(|(i, s)| &state.points[*i] != s)))
            {
                return Err(Error::Input(
                    "JOINT BIRTH frozen generation mismatch".into(),
                ));
            }
            if d.birth_policy.is_some() {
                for &i in &indices {
                    if d.birth_policy == Some(BirthPolicy::StressFree)
                        && state.points[i]
                            .birth
                            .as_ref()
                            .is_some_and(|b| b.initial_traction_kpa != [0.; 3])
                    {
                        return Err(Error::Input(
                            "mode B cannot retain an initial contact force".into(),
                        ));
                    }
                    if (l.generation > 0) != state.points[i].birth.is_some() {
                        return Err(Error::Input(
                            "JOINT BIRTH missing/unexpected birth reference".into(),
                        ));
                    }
                }
            }
        }
        Ok(())
    }
    pub fn tied_pairs(&self, state: &State, active: &[bool]) -> Result<Vec<[usize; 2]>> {
        self.validate_lifecycle(state)?;
        let mut pairs = BTreeSet::new();
        for p in &self.points {
            if Self::active(p, active) && state.lifecycle.get(p.slot(state)).is_some_and(|l| l.tied)
            {
                for (&a, &b) in p.info.nodes_a.iter().zip(&p.info.nodes_b) {
                    if a != b {
                        pairs.insert([a, b]);
                    }
                }
            }
        }
        Ok(pairs.into_iter().collect())
    }
    fn selected_joints(&self, ids: &[String]) -> Result<BTreeSet<usize>> {
        if ids.is_empty() || ids.iter().collect::<BTreeSet<_>>().len() != ids.len() {
            return Err(Error::Input(
                "JOINT BIRTH/DEATH requires unique joint ids".into(),
            ));
        }
        ids.iter()
            .map(|id| {
                self.configuration
                    .joints
                    .iter()
                    .position(|d| &d.id == id && d.birth_policy.is_some())
                    .ok_or_else(|| {
                        Error::Input(
                            "JOINT BIRTH/DEATH needs a known joint with explicit birth_policy"
                                .into(),
                        )
                    })
            })
            .collect()
    }
    fn face_region(&self, slot: usize) -> FaceRegion {
        let points = &self.points[self.faces[slot].clone()];
        let p = &points[0].info;
        FaceRegion {
            face_index: p.face_index,
            owner_a: p.owner_a,
            owner_b: p.owner_b,
            area_m2: points.iter().map(|p| p.info.weight_m2).sum(),
        }
    }

    fn frozen_work(&self, points: &[(usize, PointState)]) -> f64 {
        points
            .iter()
            .map(|(i, s)| self.points[*i].info.weight_m2 * s.endpoint_work_kpa_m)
            .sum()
    }

    // Convert the legacy whole-interface ledger only at a domain change.
    // The fixed mesh partitions its old quadrature rows without interpolating history.
    fn split_faces(&self, state: &mut State) {
        if state.lifecycle.first().is_none_or(|l| l.region.is_some()) {
            return;
        }
        state.lifecycle = self
            .faces
            .iter()
            .enumerate()
            .map(|(slot, indices)| {
                let old = &state.lifecycle[self.points[indices.start].material];
                let mut face = old.clone();
                face.region = Some(self.face_region(slot));
                for a in &mut face.archived {
                    a.points.retain(|(i, _)| indices.contains(i));
                    a.endpoint_work_kn_m = Some(self.frozen_work(&a.points));
                }
                face
            })
            .collect();
    }

    fn freeze(&self, state: &mut State, slot: usize) {
        let points: Vec<_> = self
            .points
            .iter()
            .enumerate()
            .filter(|(_, p)| p.slot(state) == slot)
            .map(|(i, _)| (i, state.points[i].clone()))
            .collect();
        let l = &mut state.lifecycle[slot];
        l.archived.push(Generation {
            generation: l.generation,
            endpoint_work_kn_m: l.region.as_ref().map(|_| self.frozen_work(&points)),
            points,
        });
        l.active = false;
        l.tied = true;
    }
    pub fn death(&self, state: &State, ids: &[String]) -> Result<State> {
        self.validate_lifecycle(state)?;
        let mut next = state.clone();
        for di in self.selected_joints(ids)? {
            let slots: BTreeSet<_> = self
                .points
                .iter()
                .filter(|p| p.material == di && state.lifecycle[p.slot(state)].active)
                .map(|p| p.slot(state))
                .collect();
            if slots.is_empty() {
                return Err(Error::Input("JOINT DEATH requires active contact".into()));
            }
            for slot in slots {
                self.freeze(&mut next, slot);
            }
        }
        Ok(next)
    }

    pub fn freeze_removed_faces(
        &self,
        state: &State,
        before: &[bool],
        after: &[bool],
    ) -> Result<State> {
        self.validate_lifecycle(state)?;
        if before.len() != self.element_count || after.len() != self.element_count {
            return Err(Error::Input("JOINT domain layout mismatch".into()));
        }
        let removed: Vec<_> = self
            .points
            .iter()
            .filter(|p| {
                self.configuration.joints[p.material].birth_policy.is_some()
                    && Self::active(p, before)
                    && !Self::active(p, after)
                    && state.lifecycle[p.slot(state)].active
            })
            .collect();
        let mut next = state.clone();
        if !removed.is_empty() {
            self.split_faces(&mut next);
            for slot in removed
                .iter()
                .map(|p| p.slot(&next))
                .collect::<BTreeSet<_>>()
            {
                self.freeze(&mut next, slot);
            }
        }
        Ok(next)
    }

    /// Domain activation supplies geometry. Explicit JOINT BIRTH supplies the policy.
    pub fn prepare_body_birth(
        &self,
        state: &State,
        before: &[bool],
        after: &[bool],
    ) -> Result<State> {
        self.validate_lifecycle(state)?;
        if before.len() != self.element_count || after.len() != self.element_count {
            return Err(Error::Input("JOINT domain layout mismatch".into()));
        }
        let mut next = state.clone();
        for (di, d) in self.configuration.joints.iter().enumerate() {
            let added: Vec<_> = self
                .points
                .iter()
                .filter(|p| p.material == di && !Self::active(p, before) && Self::active(p, after))
                .collect();
            if added.is_empty() {
                continue;
            }
            if d.birth_policy.is_none() {
                return Err(Error::Input(
                    "body BIRTH requires explicit JOINT birth_policy".into(),
                ));
            }
            let partial = self
                .points
                .iter()
                .any(|p| p.material == di && (Self::active(p, before) || !Self::active(p, after)));
            if d.water.is_some() && partial {
                return Err(Error::Input("partial JOINT BIRTH is unverified".into()));
            }
            if partial || d.crack_tip.is_some() {
                self.split_faces(&mut next);
            }
            for slot in added.iter().map(|p| p.slot(&next)).collect::<BTreeSet<_>>() {
                if next.lifecycle[slot].active {
                    self.freeze(&mut next, slot);
                }
                next.lifecycle[slot].tied = true;
            }
        }
        Ok(next)
    }
    /// The input is F_bulk + F_exterior + F_joint,water - I_material - I_pore
    /// minus every other effective interface force. Only the work-conjugate
    /// unconstrained components determine a constraint transfer multiplier.
    pub fn birth(
        &self,
        state: &State,
        u: &[[f64; 3]],
        active: &[bool],
        ids: &[String],
        residual: &[f64],
        free_rows: &[Vec<(usize, f64)>],
    ) -> Result<State> {
        self.validate_layout(state, u, active)?;
        if residual.len() != 3 * u.len() || free_rows.len() != residual.len() {
            return Err(Error::Input(
                "JOINT BIRTH force recovery layout mismatch".into(),
            ));
        }
        let selected = self.selected_joints(ids)?;
        let mut next = state.clone();
        if selected.iter().any(|&di| {
            self.configuration.joints[di].water.is_none()
                && self
                    .points
                    .iter()
                    .any(|p| p.material == di && !Self::active(p, active))
        }) {
            self.split_faces(&mut next);
        }
        let mut transfer = Vec::new();
        let mut born_points = Vec::new();
        for &di in &selected {
            let mut added = false;
            for (i, p) in self
                .points
                .iter()
                .enumerate()
                .filter(|(_, p)| p.material == di)
            {
                let life = &next.lifecycle[p.slot(&next)];
                if !Self::active(p, active) && life.region.is_some() {
                    continue;
                }
                if life.active && life.region.is_some() {
                    continue;
                }
                if life.active || !life.tied || !Self::active(p, active) {
                    return Err(Error::Input(
                        "JOINT BIRTH needs inactive contact and both active bodies".into(),
                    ));
                }
                for (&a, &b) in p.info.nodes_a.iter().zip(&p.info.nodes_b) {
                    if (Vector3::from(u[a]) - Vector3::from(u[b])).norm() > 1e-10 {
                        return Err(Error::Input(
                            "JOINT BIRTH pre-birth equal-displacement constraint is not satisfied"
                                .into(),
                        ));
                    }
                }
                added = true;
                born_points.push(i);
                if life.policy.is_some_and(BirthPolicy::transfers_force) {
                    transfer.push(i);
                }
            }
            if !added {
                return Err(Error::Input(
                    "JOINT BIRTH needs inactive contact and both active bodies".into(),
                ));
            }
        }
        let mut traction = BTreeMap::new();
        if !transfer.is_empty() {
            // Recover the nodal tie multipliers BEFORE fitting a GP field.
            // Boundary reactions are eliminated by the untied nullspace;
            // redundant prescribed components use the minimum-norm multiplier.
            // The second solve then checks ALL physical nodal components.
            let nodes: BTreeSet<_> = transfer
                .iter()
                // Keep shared tip DOFs as zero contact-force rows. Removing
                // their zero coefficients from connectivity must not remove
                // their free residual from the exact recovery/rejection gate.
                .flat_map(|&i| {
                    self.points[i]
                        .info
                        .nodes_a
                        .iter()
                        .chain(&self.points[i].info.nodes_b)
                        .copied()
                })
                .collect();
            let physical: BTreeMap<_, _> = nodes
                .iter()
                .flat_map(|&n| (0..3).map(move |a| 3 * n + a))
                .enumerate()
                .map(|(i, g)| (g, i))
                .collect();
            let free: BTreeSet<_> = physical
                .keys()
                .flat_map(|&g| free_rows[g].iter().map(|&(r, _)| r))
                .collect();
            let free_index: BTreeMap<_, _> =
                free.iter().enumerate().map(|(i, &r)| (r, i)).collect();
            let pairs: BTreeSet<_> = transfer
                .iter()
                .flat_map(|&i| {
                    self.points[i]
                        .info
                        .nodes_a
                        .iter()
                        .zip(&self.points[i].info.nodes_b)
                        // A shared tip is one physical node, not a tie
                        // multiplier. Its jump and contact force cancel exactly.
                        .filter(|(a, b)| a != b)
                        .map(|(&a, &b)| [a, b])
                })
                .collect();
            if transfer.len() > 2000 || physical.len() > 6000 || pairs.len() > 2000 {
                return Err(Error::Input(
                    "mode A exact projection budget exceeded".into(),
                ));
            }
            if free.is_empty() {
                return Err(Error::Input(
                    "mode A transfer is indeterminate under prescribed interface displacements"
                        .into(),
                ));
            }
            let mut pair_matrix = DMatrix::<f64>::zeros(free.len(), 3 * pairs.len());
            let mut free_rhs = DVector::<f64>::zeros(free.len());
            let mut outside = BTreeMap::<usize, f64>::new();
            // T3 checks the complement of the selected patch too: contact must
            // supply zero there. Legacy whole-interface cases retain their
            // interface projection and the existing bulk equilibrium criterion.
            let certify_global = transfer.iter().any(|&i| {
                let p = &self.points[i];
                self.configuration.joints[p.material].crack_tip.is_some()
                    || next.lifecycle[p.slot(&next)].region.is_some()
            });
            for (g, terms) in free_rows.iter().enumerate() {
                if !certify_global && !physical.contains_key(&g) {
                    continue;
                }
                for &(r, w) in terms {
                    if let Some(&row) = free_index.get(&r) {
                        free_rhs[row] += w * residual[g];
                    } else {
                        *outside.entry(r).or_default() += w * residual[g];
                    }
                }
            }
            if let Some((row, force)) = outside
                .iter()
                .find(|(_, v)| !v.is_finite() || v.abs() > 1e-10)
            {
                return Err(Error::Input(format!(
                    "mode A constraint transfer does not reproduce every free nodal residual; free row={row} requested_kn={force} reproduced_kn=0"
                )));
            }
            for (j, &[a, b]) in pairs.iter().enumerate() {
                for (node, sign) in [(a, -1.), (b, 1.)] {
                    for axis in 0..3 {
                        for &(r, w) in &free_rows[3 * node + axis] {
                            pair_matrix[(free_index[&r], 3 * j + axis)] += sign * w;
                        }
                    }
                }
            }
            let multipliers = recover_forces(
                &pair_matrix,
                &free_rhs,
                &free.iter().copied().collect::<Vec<_>>(),
                "mode A constraint transfer does not reproduce every free nodal residual",
            )?;
            let mut rhs = DVector::<f64>::zeros(physical.len());
            for (j, &[a, b]) in pairs.iter().enumerate() {
                for axis in 0..3 {
                    rhs[physical[&(3 * a + axis)]] -= multipliers[3 * j + axis];
                    rhs[physical[&(3 * b + axis)]] += multipliers[3 * j + axis];
                }
            }
            let mut matrix = DMatrix::<f64>::zeros(physical.len(), 3 * transfer.len());
            for (j, &pi) in transfer.iter().enumerate() {
                let p = &self.points[pi];
                for &(node, c) in &p.coefficients {
                    for a in 0..3 {
                        matrix[(physical[&(3 * node + a)], 3 * j + a)] +=
                            p.info.weight_m2.sqrt() * c;
                    }
                }
            }
            let rows: Vec<_> = physical.keys().copied().collect();
            let failure =
                "mode A GP initial forces cannot reproduce every physical nodal transfer force";
            let mean_faces: BTreeSet<_> = transfer
                .iter()
                .filter_map(|&i| {
                    let p = &self.points[i];
                    (self.configuration.joints[p.material].birth_policy
                        == Some(BirthPolicy::TransferForceFaceMean))
                    .then_some((p.material, p.info.face_index))
                })
                .collect();
            let q = if mean_faces.is_empty() {
                recover_forces(&matrix, &rhs, &rows, failure)?
            } else {
                // Lift a constant traction per new face to the weighted GP space.
                // Fit its virtual work on the physical nodes, then recover the exact
                // remaining force in the same space as the legacy transfer.
                let columns: BTreeMap<_, _> = mean_faces
                    .into_iter()
                    .enumerate()
                    .map(|(i, k)| (k, i))
                    .collect();
                let mut lift = DMatrix::zeros(matrix.ncols(), 3 * columns.len());
                for (j, &i) in transfer.iter().enumerate() {
                    let p = &self.points[i];
                    if let Some(&f) = columns.get(&(p.material, p.info.face_index)) {
                        for axis in 0..3 {
                            lift[(3 * j + axis, 3 * f + axis)] = p.info.weight_m2.sqrt();
                        }
                    }
                }
                let mean = (&matrix * &lift)
                    .svd(true, true)
                    .solve(&rhs, 1e-12)
                    .map_err(|err| Error::Input(format!("face mean transfer: {err}")))?;
                let seed = lift * mean;
                let correction =
                    recover_forces(&matrix, &(&rhs - &matrix * &seed), &rows, failure)?;
                let field = seed + correction;
                let error = &matrix * &field - &rhs;
                if error.iter().any(|v| !v.is_finite() || v.abs() > 1e-10) {
                    return Err(Error::Input(failure.into()));
                }
                field
            };
            let reproduced = &matrix * &q;
            let mut violations = Vec::new();
            for (j, &pi) in transfer.iter().enumerate() {
                let scale = self.points[pi].info.weight_m2.sqrt();
                let value = [q[3 * j] / scale, q[3 * j + 1] / scale, q[3 * j + 2] / scale];
                let point = &self.points[pi];
                let material = &self.configuration.joints[point.material].material;
                if let Err(err) =
                    material.birth(point.info.normal, Self::jump(point, u).into(), value)
                {
                    violations.push(format!(
                        "GP={pi} face={} q_kpa={value:?} {err}",
                        point.info.face_index
                    ));
                }
                traction.insert(pi, value);
            }
            if !violations.is_empty() {
                return Err(Error::Input(format!(
                    "JOINT BIRTH rejected trial: {}; physical_dofs={:?} requested_kn={:?} reproduced_kn={:?}",
                    violations.join(", "),
                    physical.keys().collect::<Vec<_>>(),
                    rhs.as_slice(),
                    reproduced.as_slice()
                )));
            }
        }
        let mut slots = BTreeSet::new();
        for i in born_points {
            let p = &self.points[i];
            slots.insert(p.slot(&next));
            let q = traction.get(&i).copied().unwrap_or([0.; 3]);
            let material = &self.configuration.joints[p.material].material;
            next.points[i] = material
                .birth(p.info.normal, Self::jump(p, u).into(), q)
                .map_err(|err| {
                    Error::Input(format!(
                        "{err}; trial GP={i} face={} q_kpa={q:?}",
                        p.info.face_index
                    ))
                })?;
        }
        for slot in slots {
            let l = &mut next.lifecycle[slot];
            l.active = true;
            l.tied = false;
            l.generation = l
                .generation
                .checked_add(1)
                .ok_or_else(|| Error::Input("JOINT generation overflow".into()))?;
        }
        self.validate_state(&next, u, active)?;
        Ok(next)
    }
}

fn recover_forces(
    matrix: &DMatrix<f64>,
    rhs: &DVector<f64>,
    rows: &[usize],
    failure: &str,
) -> Result<DVector<f64>> {
    const CUTOFF: f64 = 1e-12;
    let svd = matrix.clone().svd(true, true);
    let solution = svd
        .solve(rhs, CUTOFF)
        .map_err(|e| Error::Input(format!("{failure}: {e}")))?;
    let reproduced = matrix * &solution;
    if reproduced
        .iter()
        .zip(rhs)
        .any(|(a, b)| !(a - b).is_finite() || (a - b).abs() > 1e-10)
    {
        return Err(Error::Input(format!(
            "{failure}; rows={rows:?} requested_kn={:?} reproduced_kn={:?} rank={} singular_values={:?} cutoff={CUTOFF}",
            rhs.as_slice(),
            reproduced.as_slice(),
            svd.rank(CUTOFF),
            svd.singular_values.as_slice()
        )));
    }
    Ok(solution)
}
