//! Opt-in, budgeted exact condensation. The legacy500-node path stays unchanged.
use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StorageBudget {
    pub maximum_nodes: usize,
    pub maximum_compilation_mib: usize,
    pub maximum_cached_mib: usize,
}

const MIB: usize = 1024 * 1024;
type SparseColumn = Vec<(usize, f64)>;

impl StorageBudget {
    pub(super) fn validate(&self) -> Result<()> {
        if !(1..=1000).contains(&self.maximum_nodes)
            || !(1..=2048).contains(&self.maximum_compilation_mib)
            || self.maximum_cached_mib > 1024
        {
            return Err(Error::Input(
                "invalid explicit structural storage budget".into(),
            ));
        }
        Ok(())
    }
}

fn storage_estimate(nd: usize, ni: usize, nf: usize) -> Result<(usize, usize)> {
    let overflow = || Error::Input("structural matrix storage estimate overflow".into());
    let product = |a: usize, b: usize| a.checked_mul(b).ok_or_else(overflow);
    let sum = |items: &[usize]| {
        items
            .iter()
            .try_fold(0usize, |a, b| a.checked_add(*b).ok_or_else(overflow))
    };
    let kk = product(nd, nd)?;
    let tt = product(nd, ni)?;
    let qq = product(nd, nf)?;
    let ff = product(nf, nf)?;
    let ii = product(ni, ni)?;
    // K, T, Q, Kqq, its Cholesky storage, Kc and the kinematic offset.
    let resident = sum(&[kk, tt, qq, ff, ff, ii, nd])?;
    // Conservative extra buffers for Kqt, solve, Schur products and factorization.
    // Existing cache storage and other FEM/model data are separate budgets.
    let fi = product(nf, ni)?;
    let transient = sum(&[fi, fi, fi, ii, ii, ii, ff])?;
    Ok((
        product(resident, 8)?,
        product(sum(&[resident, transient])?, 8)?,
    ))
}

fn project(k: &DMatrix<f64>, left: &[SparseColumn], right: &[SparseColumn]) -> DMatrix<f64> {
    DMatrix::from_fn(left.len(), right.len(), |i, j| {
        right[j]
            .iter()
            .map(|&(b, vb)| left[i].iter().map(|&(a, va)| va * k[(a, b)]).sum::<f64>() * vb)
            .sum()
    })
}

impl Model {
    pub(super) fn compiled_bounded(&self, state: &State, key: Vec<bool>) -> Result<Rc<Compiled>> {
        let budget = self.configuration.storage_budget.as_ref().unwrap();
        let n = self.positions.len();
        let nd = 6 * n;
        let mut active_nodes = vec![false; n];
        for (i, m) in self.members.iter().enumerate() {
            if state.active_members[i] {
                for &node in &m.nodes {
                    active_nodes[node] = true;
                }
            }
        }
        for (i, p) in self.plates.iter().enumerate() {
            if state.active_plates[i] {
                for &node in &p.nodes {
                    active_nodes[node] = true;
                }
            }
        }
        let mesh_nodes: BTreeSet<_> = self
            .attachments
            .iter()
            .enumerate()
            .filter(|(i, _)| active_nodes[*i])
            .filter_map(|(_, a)| *a)
            .collect();
        let interface_dofs: Vec<_> = mesh_nodes
            .iter()
            .flat_map(|&v| [3 * v, 3 * v + 1, 3 * v + 2])
            .collect();
        let interface_ids: BTreeMap<_, _> = interface_dofs
            .iter()
            .copied()
            .enumerate()
            .map(|(i, v)| (v, i))
            .collect();
        let mut t_columns: Vec<SparseColumn> = vec![Vec::new(); interface_dofs.len()];
        let mut columns: Vec<SparseColumn> = Vec::new();
        let mut offset = DVector::zeros(nd);
        let mut normals = Vec::new();
        for i in 0..n {
            if !active_nodes[i] {
                continue;
            }
            for part in 0..2 {
                if part == 0
                    && let Some(node) = self.attachments[i]
                {
                    for a in 0..3 {
                        t_columns[interface_ids[&(3 * node + a)]].push((6 * i + a, 1.));
                    }
                    continue;
                }
                let constraints = if part == 0 {
                    &self.configuration.translation_constraints
                } else {
                    &self.configuration.rotation_constraints
                };
                let mut equations = Vec::new();
                for c in constraints
                    .iter()
                    .filter(|c| c.node_id == self.configuration.nodes[i].id)
                {
                    let d = Vector3::from(c.direction);
                    equations.push((d / d.norm(), c.value / d.norm()));
                }
                let (u0, free, fixed) = basis(&equations)?;
                for a in 0..3 {
                    offset[6 * i + 3 * part + a] = u0[a];
                }
                for normal in fixed {
                    normals.push((6 * i + 3 * part, normal));
                }
                for b in free {
                    let scale = if part == 1 {
                        self.configuration.rotation_length_scale_m
                    } else {
                        1.
                    };
                    columns.push(
                        (0..3)
                            .filter_map(|a| {
                                let value = b[a] / scale;
                                (value != 0.).then_some((6 * i + 3 * part + a, value))
                            })
                            .collect(),
                    );
                }
            }
        }
        if columns.len() > self.configuration.maximum_condensed_dofs {
            return Err(Error::Input("beam interior DOF budget exceeded".into()));
        }
        let (resident, peak) = storage_estimate(nd, t_columns.len(), columns.len())?;
        if peak > budget.maximum_compilation_mib * MIB {
            return Err(Error::Input(format!(
                "structural compilation matrix estimate {peak} bytes exceeds explicit storage budget before allocation"
            )));
        }
        let cache_limit = budget.maximum_cached_mib * MIB;
        {
            let mut cache = self.cache.borrow_mut();
            let cached = cache
                .values()
                .try_fold(0usize, |a, c| a.checked_add(c.matrix_storage_bytes))
                .ok_or_else(|| Error::Input("structural cache storage overflow".into()))?;
            if cached.checked_add(resident).is_none_or(|v| v > cache_limit) {
                cache.clear();
            }
        }
        // No dense quadratic matrix is allocated before the estimate above.
        let mut k = DMatrix::zeros(nd, nd);
        for (i, m) in self.members.iter().enumerate() {
            if state.active_members[i] {
                for a in 0..12 {
                    for b in 0..12 {
                        k[(6 * m.nodes[a / 6] + a % 6, 6 * m.nodes[b / 6] + b % 6)] +=
                            m.kernel.stiffness[(a, b)];
                    }
                }
            }
        }
        for (i, p) in self.plates.iter().enumerate() {
            if state.active_plates[i] {
                for a in 0..18 {
                    for b in 0..18 {
                        k[(6 * p.nodes[a / 6] + a % 6, 6 * p.nodes[b / 6] + b % 6)] +=
                            p.kernel.stiffness[(a, b)];
                    }
                }
            }
        }
        let mut q = DMatrix::zeros(nd, columns.len());
        for (j, column) in columns.iter().enumerate() {
            for &(i, v) in column {
                q[(i, j)] = v;
            }
        }
        let mut t = DMatrix::zeros(nd, t_columns.len());
        for (j, column) in t_columns.iter().enumerate() {
            for &(i, v) in column {
                t[(i, j)] = v;
            }
        }
        let kqq = project(&k, &columns, &columns);
        let kqt = project(&k, &columns, &t_columns);
        let chol = if columns.is_empty() {
            None
        } else {
            let c = kqq.clone().cholesky().ok_or_else(|| {
                Error::Solve(
                    "unconstrained beam mechanism; inspect torsion/releases/supports".into(),
                )
            })?;
            let maximum = (0..kqq.nrows())
                .map(|i| kqq[(i, i)].abs())
                .fold(0., f64::max);
            if c.l().diagonal().iter().any(|d| d * d <= 1e-12 * maximum) {
                return Err(Error::Solve(
                    "singular/ill-conditioned beam interior; no stabilization applied".into(),
                ));
            }
            Some(c)
        };
        let kc = project(&k, &t_columns, &t_columns)
            - chol.as_ref().map_or_else(
                || DMatrix::zeros(t.ncols(), t.ncols()),
                |c| kqt.transpose() * c.solve(&kqt),
            );
        let actual = (k.len()
            + t.len()
            + q.len()
            + kqq.len()
            + kc.len()
            + offset.len()
            + if chol.is_some() {
                columns.len() * columns.len()
            } else {
                0
            })
            * 8;
        if actual != resident {
            return Err(Error::Input(
                "structural resident matrix accounting mismatch".into(),
            ));
        }
        let compiled = Rc::new(Compiled {
            active_nodes,
            k,
            t,
            q,
            offset,
            chol,
            kqq,
            kc,
            interface_dofs,
            normals,
            matrix_storage_bytes: actual,
        });
        if actual <= cache_limit {
            self.cache.borrow_mut().insert(key, compiled.clone());
        }
        Ok(compiled)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn cache_model(cache: usize) -> Result<Model> {
        let mut mesh: Mesh =
            serde_json::from_str(include_str!("../../examples/mc-column-mesh.json"))?;
        mesh.node_sets.insert("Base".into(), vec![0]);
        let mut nodes = vec![
            serde_json::json!({"id":"base","position_m":mesh.nodes[0],"soil_node_set":"Base"}),
            serde_json::json!({"id":"tip","position_m":[2.,0.,0.],"soil_node_set":null}),
        ];
        for i in 2..40 {
            nodes.push(serde_json::json!({"id":format!("unused{i}"),"position_m":mesh.nodes[0],"soil_node_set":"Base"}));
        }
        let member = |id: &str| serde_json::json!({"id":id,"nodes":["base","tip"],"section_id":"steel","local_y":[0.,1.,0.],"initially_active":true});
        let cfg: Configuration = serde_json::from_value(
            serde_json::json!({"schema":"soil-fem-beams/1","nodes":nodes,"sections":[{"id":"steel","young_kpa":2e8,"poisson":0.3,"area_m2":0.02,"iy_m4":0.0001,"iz_m4":0.0002,"torsion_j_m4":0.00005,"unit_weight_kn_m3":0.,"formulation":{"kind":"euler_bernoulli"}}],"members":[member("a"),member("b")],"rotation_length_scale_m":2.,"maximum_condensed_dofs":100,"storage_budget":{"maximum_nodes":500,"maximum_compilation_mib":64,"maximum_cached_mib":cache},"rotation_constraints":[{"node_id":"base","direction":[1.,0.,0.],"value":0.},{"node_id":"base","direction":[0.,1.,0.],"value":0.},{"node_id":"base","direction":[0.,0.,1.],"value":0.}]}),
        )?;
        Model::new(&mesh, cfg)
    }
    #[test]
    fn uncached_and_evicted_compilations_reproduce_identical_schur_matrices() -> Result<()> {
        let disabled = cache_model(0)?;
        let state = disabled.initial_state();
        let a = disabled.compiled(&state)?;
        let b = disabled.compiled(&state)?;
        assert!(!Rc::ptr_eq(&a, &b));
        assert_eq!(a.kc, b.kc);
        assert!(disabled.cache.borrow().is_empty());
        let model = cache_model(1)?;
        let mut state = model.initial_state();
        state.active_members = vec![true, false];
        let first = model.compiled(&state)?;
        let expected = first.kc.clone();
        drop(first);
        state.active_members = vec![false, true];
        model.compiled(&state)?;
        state.active_members = vec![true, true];
        model.compiled(&state)?;
        assert_eq!(model.cache.borrow().len(), 1);
        assert!(
            model
                .cache
                .borrow()
                .values()
                .map(|c| c.matrix_storage_bytes)
                .sum::<usize>()
                <= MIB
        );
        state.active_members = vec![true, false];
        assert_eq!(model.compiled(&state)?.kc, expected);
        Ok(())
    }
    #[test]
    fn sparse_constraint_projection_matches_independent_dense_products() {
        let a = DMatrix::from_fn(12, 12, |i, j| ((i * 17 + j * 7 + 3) % 23) as f64 / 23.);
        let k = a.transpose() * a + DMatrix::identity(12, 12);
        let q = vec![
            vec![(3, 0.6), (4, -0.8)],
            vec![(4, 0.8), (3, 0.6)],
            vec![(9, 0.25), (10, -0.75), (11, 0.5)],
        ];
        let t = vec![
            vec![(0, 1.), (6, 1.)],
            vec![(1, 1.)],
            vec![(2, 1.), (8, 1.)],
        ];
        let dense = |cols: &Vec<SparseColumn>| {
            DMatrix::from_fn(12, cols.len(), |i, j| {
                cols[j].iter().find(|x| x.0 == i).map_or(0., |x| x.1)
            })
        };
        let (qd, td) = (dense(&q), dense(&t));
        for (got, want) in [
            (project(&k, &q, &q), qd.transpose() * &k * &qd),
            (project(&k, &q, &t), qd.transpose() * &k * &td),
            (project(&k, &t, &t), td.transpose() * &k * &td),
        ] {
            assert!((&got - &want).norm() < 1e-13 * want.norm());
        }
    }
    #[test]
    fn storage_arithmetic_overflow_and_invalid_budgets_are_rejected() {
        assert!(storage_estimate(usize::MAX, 3, 3).is_err());
        for (nodes, compilation, cache) in [
            (1001, 1024, 512),
            (0, 1024, 512),
            (750, 0, 512),
            (750, 2049, 512),
            (750, 1024, 1025),
        ] {
            assert!(
                StorageBudget {
                    maximum_nodes: nodes,
                    maximum_compilation_mib: compilation,
                    maximum_cached_mib: cache
                }
                .validate()
                .is_err()
            );
        }
        let (resident, peak) = storage_estimate(3972, 1986, 1986).unwrap();
        assert!(
            resident > 300 * MIB && resident < 512 * MIB && peak > resident && peak < 1024 * MIB
        );
    }
}
