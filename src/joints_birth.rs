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
        if state.lifecycle.len() != self.configuration.joints.len() {
            return Err(Error::Input("JOINT BIRTH lifecycle layout mismatch".into()));
        }
        for (di, (d, l)) in self
            .configuration
            .joints
            .iter()
            .zip(&state.lifecycle)
            .enumerate()
        {
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
                .filter(|(_, p)| p.material == di)
                .map(|(i, _)| i)
                .collect();
            let mut previous = 0;
            for a in &l.archived {
                if a.generation == 0
                    || a.generation <= previous
                    || a.generation > l.generation
                    || a.points.iter().map(|(i, _)| *i).collect::<Vec<_>>() != indices
                {
                    return Err(Error::Input(
                        "JOINT BIRTH archived generation mismatch".into(),
                    ));
                }
                for (_, s) in &a.points {
                    d.material.validate_state(s)?;
                }
                previous = a.generation;
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
            if Self::active(p, active) && state.lifecycle.get(p.material).is_some_and(|l| l.tied) {
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
    fn archive(&self, state: &mut State, di: usize) {
        let l = &mut state.lifecycle[di];
        if l.generation > 0
            && l.archived
                .last()
                .is_none_or(|a| a.generation != l.generation)
        {
            l.archived.push(Generation {
                generation: l.generation,
                points: self
                    .points
                    .iter()
                    .enumerate()
                    .filter(|(_, p)| p.material == di)
                    .map(|(i, _)| (i, state.points[i].clone()))
                    .collect(),
            });
        }
    }
    pub fn death(&self, state: &State, ids: &[String]) -> Result<State> {
        self.validate_lifecycle(state)?;
        let mut next = state.clone();
        for di in self.selected_joints(ids)? {
            if !next.lifecycle[di].active {
                return Err(Error::Input("JOINT DEATH requires active contact".into()));
            }
            self.archive(&mut next, di);
            next.lifecycle[di].active = false;
            next.lifecycle[di].tied = true;
        }
        Ok(next)
    }
    /// After a dry body BIRTH the interface remains tied until its own BIRTH.
    pub fn prepare_body_birth(
        &self,
        state: &State,
        before: &[bool],
        after: &[bool],
    ) -> Result<State> {
        self.validate_lifecycle(state)?;
        let mut next = state.clone();
        for di in 0..self.configuration.joints.len() {
            let changes = self
                .points
                .iter()
                .any(|p| p.material == di && !Self::active(p, before) && Self::active(p, after));
            if changes {
                if self.configuration.joints[di].birth_policy.is_none() {
                    return Err(Error::Input(
                        "body BIRTH requires explicit JOINT birth_policy".into(),
                    ));
                }
                if self
                    .points
                    .iter()
                    .any(|p| p.material == di && !Self::active(p, after))
                {
                    return Err(Error::Input("partial JOINT BIRTH is unverified".into()));
                }
                self.archive(&mut next, di);
                next.lifecycle[di].active = false;
                next.lifecycle[di].tied = true;
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
        let mut transfer = Vec::new();
        for &di in &selected {
            let l = &state.lifecycle[di];
            if l.active
                || !l.tied
                || self
                    .points
                    .iter()
                    .any(|p| p.material == di && !Self::active(p, active))
            {
                return Err(Error::Input(
                    "JOINT BIRTH needs inactive contact and both active bodies".into(),
                ));
            }
            for p in self.points.iter().filter(|p| p.material == di) {
                for (&a, &b) in p.info.nodes_a.iter().zip(&p.info.nodes_b) {
                    if (Vector3::from(u[a]) - Vector3::from(u[b])).norm() > 1e-10 {
                        return Err(Error::Input(
                            "JOINT BIRTH pre-birth equal-displacement constraint is not satisfied"
                                .into(),
                        ));
                    }
                }
            }
            if l.policy == Some(BirthPolicy::TransferForce) {
                if !l.tied {
                    return Err(Error::Input(
                        "mode A requires a pre-birth equal-displacement constraint".into(),
                    ));
                }
                transfer.extend(
                    self.points
                        .iter()
                        .enumerate()
                        .filter(|(_, p)| p.material == di)
                        .map(|(i, _)| i),
                );
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
                .flat_map(|&i| self.points[i].coefficients.iter().map(|&(n, _)| n))
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
            for &g in physical.keys() {
                for &(r, w) in &free_rows[g] {
                    free_rhs[free_index[&r]] += w * residual[g];
                }
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
            let multipliers = pair_matrix
                .clone()
                .svd(true, true)
                .solve(&free_rhs, 1e-12)
                .map_err(|_| Error::Input("mode A constraint transfer recovery failed".into()))?;
            if (&pair_matrix * &multipliers - free_rhs)
                .iter()
                .any(|v| !v.is_finite() || v.abs() > 1e-10)
            {
                return Err(Error::Input(
                    "mode A constraint transfer does not reproduce every free nodal residual"
                        .into(),
                ));
            }
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
            let q = matrix
                .clone()
                .svd(true, true)
                .solve(&rhs, 1e-12)
                .map_err(|_| Error::Input("mode A nodal transfer projection failed".into()))?;
            if (&matrix * &q - &rhs)
                .iter()
                .any(|v| !v.is_finite() || v.abs() > 1e-10)
            {
                return Err(Error::Input(
                    "mode A GP initial forces cannot reproduce every physical nodal transfer force"
                        .into(),
                ));
            }
            for (j, &pi) in transfer.iter().enumerate() {
                let scale = self.points[pi].info.weight_m2.sqrt();
                traction.insert(
                    pi,
                    [q[3 * j] / scale, q[3 * j + 1] / scale, q[3 * j + 2] / scale],
                );
            }
        }
        for di in selected {
            self.archive(&mut next, di);
            for (i, p) in self
                .points
                .iter()
                .enumerate()
                .filter(|(_, p)| p.material == di)
            {
                next.points[i] = self.configuration.joints[di].material.birth(
                    p.info.normal,
                    Self::jump(p, u).into(),
                    traction.get(&i).copied().unwrap_or([0.; 3]),
                )?;
            }
            let l = &mut next.lifecycle[di];
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
