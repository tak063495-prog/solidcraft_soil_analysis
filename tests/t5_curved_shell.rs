use soil_fem_3d::{
    Result,
    curved_shell::{
        Kernel, Q36, Tying,
        geometry::{self, Face, Geometry, V},
    },
    plates::{Section, offset},
};
fn section() -> Section {
    Section {
        id: "t5".into(),
        young_kpa: 20000.,
        poisson: 0.25,
        thickness_m: 0.1,
        unit_weight_kn_m3: 0.,
        drilling_factor: 1e-4,
    }
}
fn face() -> Face {
    Face {
        id: "f".into(),
        group: "smooth".into(),
        nodes: [0, 1, 2, 3, 4, 5],
        offset_m: 0.,
    }
}
fn patch(a: f64, b: f64, e: f64) -> Geometry {
    let x = geometry::NODES.map(|[r, s]| V::new(r, s, a * r * r + b * s * s));
    let directors = geometry::NODES.map(|[r, s]| V::new(-2. * a * r, -2. * b * s, 1.).normalize());
    Geometry { x, directors, e }
}
fn dofs(g: &Geometry, f: impl Fn(V) -> (V, V)) -> Q36 {
    Q36::from_fn(|i, _| {
        let (u, r) = f(g.x[i / 6]);
        if i % 6 < 3 { u[i % 6] } else { r[i % 6 - 3] }
    })
}
#[test]
fn curved_thickness_map_and_normalized_linearisation() -> Result<()> {
    let g = patch(0.3, -0.2, -0.05);
    let p = [0.2, 0.38];
    let s = g.surface(p)?;
    let rot = V::new(-0.2, 0.6, 0.9);
    let dd = s.variation(3, rot.cross(&g.directors[3]));
    let eps = 1e-6;
    let mut plus = g.clone();
    let mut minus = g.clone();
    plus.directors[3] += eps * rot.cross(&g.directors[3]);
    minus.directors[3] -= eps * rot.cross(&g.directors[3]);
    let (hi, lo) = (plus.surface(p)?, minus.surface(p)?);
    assert!(((hi.d - lo.d) / (2. * eps) - dd.0).norm() < 1e-9);
    for k in 0..2 {
        assert!(((hi.dd[k] - lo.dd[k]) / (2. * eps) - dd.1[k]).norm() < 1e-9);
    }
    assert!((g.basis(&s, -0.05) - g.basis(&s, 0.05)).norm() > 0.01);
    g.certify(0.1)?;
    let mut folded = g.clone();
    folded.x[3] = g.x[0] + (g.x[1] - g.x[0]) * 0.05;
    assert!(folded.certify(0.1).is_err());
    Ok(())
}
#[test]
fn director_ownership_shared_edge_and_crease() -> Result<()> {
    let mut nodes = patch(0., 0., 0.).x.map(Into::into).to_vec();
    nodes.extend([[1., 1., 0.4], [1., 0.5, 0.2], [0.5, 1., 0.2]]);
    let f = face();
    let mut h = Face {
        id: "h".into(),
        nodes: [2, 1, 6, 4, 7, 8],
        ..f.clone()
    };
    let pair = geometry::geometries(&nodes, &[f.clone(), h.clone()], &[])?;
    for (i, j) in [(1, 1), (2, 0), (4, 3)] {
        assert_eq!(pair[0].directors[i], pair[1].directors[j]);
    }
    let rev = geometry::geometries(&nodes, &[h.clone(), f.clone()], &[])?;
    assert_eq!(pair[0].directors, rev[1].directors);
    h.group = "crease".into();
    let crease = geometry::geometries(&nodes, &[f.clone(), h.clone()], &[])?;
    assert!((crease[0].directors[4] - crease[1].directors[3]).norm() > 0.1);
    h.offset_m = 0.05;
    assert!(geometry::geometries(&nodes, &[f, h], &[]).is_err());
    Ok(())
}
#[test]
fn rigid_fields_and_drilling_coefficient_series() -> Result<()> {
    let g = patch(0.3, -0.2, 0.03);
    let mut prior = None;
    for alpha in [1e-5, 1e-4, 1e-3] {
        let s = Section {
            drilling_factor: alpha,
            ..section()
        };
        let k = Kernel::new(&g, &s, Tying::default())?;
        for a in 0..3 {
            let axis = V::ith(a, 1.);
            for q in [
                dofs(&g, |_| (axis, V::zeros())),
                dofs(&g, |x| (axis.cross(&x), axis)),
            ] {
                assert!((k.stiffness * q).norm() < 1e-10);
                assert!(k.fields(&q).iter().all(|f| f.drill_rad.abs() < 1e-12));
            }
        }
        let q = Q36::from_fn(|i, _| {
            if i % 6 >= 3 {
                g.directors[i / 6][i % 6 - 3]
            } else {
                0.
            }
        });
        let work = q.dot(&(k.drill * q));
        assert!(work > 0.);
        assert!(((k.stiffness - k.drill) * q).norm() < 1e-10);
        if let Some(w) = prior {
            assert!(f64::abs(work / w - 10.) < 1e-10);
        }
        prior = Some(work);
        let eig = k.stiffness.symmetric_eigenvalues();
        assert_eq!(eig.iter().filter(|&&v| v.abs() < 1e-9).count(), 6);
    }
    Ok(())
}
#[test]
fn membrane_and_shear_tying_are_independent() -> Result<()> {
    let g = patch(0.4, 0.3, 0.);
    let config = |m, s| {
        Kernel::new(
            &g,
            &section(),
            Tying {
                membrane: m,
                shear: s,
            },
        )
    };
    let full = config(true, true)?;
    for k in [config(false, true)?, config(true, false)?] {
        assert!((k.stiffness - full.stiffness).norm() > 1e-4);
    }
    let flat = patch(0., 0., 0.);
    let k = Kernel::new(&flat, &section(), Tying::default())?;
    let q = dofs(&flat, |p| {
        (V::new(0.002 * p.x, -0.003 * p.y, 0.), V::zeros())
    });
    for f in k.fields(&q) {
        assert!((f.strain[0] - 0.002).abs() < 1e-14);
        assert!((f.strain[1] + 0.003).abs() < 1e-14);
    }
    Ok(())
}
#[test]
fn signed_offset_load_has_t4_flat_limit() -> Result<()> {
    let traction = V::new(2., -3., 1.);
    for e in [-0.05, 0., 0.05] {
        let g = patch(0., 0., e);
        let k = Kernel::new(&g, &section(), Tying::default())?;
        let flat = offset::Kernel::new(g.x.map(Into::into), &section(), e)?;
        let f = k.traction.transpose() * traction;
        let expected = flat.uniform_load(traction.into());
        assert!(
            (f.as_slice()
                .iter()
                .zip(expected.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max))
                < 1e-13
        );
        let q = dofs(&g, |x| (V::new(0.001 * x.x, 0., 0.), V::zeros()));
        let d = nalgebra::DVector::from_column_slice(q.as_slice());
        assert!((q.dot(&(k.stiffness * q)) - d.dot(&(&flat.stiffness * &d))).abs() < 1e-13);
    }
    Ok(())
}

#[test]
fn restart_replays_generations_and_all_fields_in_another_process() -> Result<()> {
    use soil_fem_3d::curved_shell::path::Case;
    use soil_fem_3d::offset_composite::path::Step;
    use std::{fs, process::Command};
    let base = patch(0.22, -0.15, 0.);
    let n = 36;
    let load = |u: f64| Step::Load {
        prescribed_configuration: (0..n)
            .map(|i| if i % 6 == 2 { u * (i / 6) as f64 } else { 0. })
            .collect(),
        free_dofs: vec![],
        nodal_load: vec![0.; n],
        plate_traction_kpa: [0.; 3],
    };
    let install = || Step::Install {
        contacts: Default::default(),
        holder: "support".into(),
    };
    let case = Case {
        soil: None,
        face_loads: Default::default(),
        nodes: base.x.map(Into::into).to_vec(),
        faces: vec![face()],
        directors: vec![],
        section: section(),
        steps: vec![
            load(1e-4),
            install(),
            load(2e-4),
            Step::Remove {
                holder: "support".into(),
            },
            install(),
            load(3e-4),
        ],
    };
    let dir = std::env::temp_dir().join(format!(
        "soil-t5-restart-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&dir)?;
    fs::write(dir.join("case"), serde_json::to_vec(&case)?)?;
    let cli = |args: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_soil-fem-3d"))
            .current_dir(&dir)
            .args(["curved-shell", "case"])
            .args(args)
            .output()
            .unwrap()
    };
    let full = cli(&["full"]);
    assert!(
        full.status.success(),
        "{}",
        String::from_utf8_lossy(&full.stderr)
    );
    let bytes = fs::read(dir.join("full"))?;
    let run: serde_json::Value = serde_json::from_slice(&bytes)?;
    assert_eq!(
        run["frames"][2]["active"]["elements"][0]["physical_energy_kn_m"],
        0.
    );
    assert_eq!(run["frames"][3]["active"], run["frames"][4]["archived"][0]);
    assert_eq!(run["frames"][5]["active"]["generation"], 2);
    for cut in [2, 3, 4, 5] {
        let cp = format!("cut-{cut}");
        assert!(cli(&[&cp, "--until", &cut.to_string()]).status.success());
        let final_file = format!("end-{cut}");
        assert!(cli(&[&final_file, "--resume", &cp]).status.success());
        assert_eq!(fs::read(dir.join(final_file))?, bytes);
    }
    let saved: serde_json::Value = serde_json::from_slice(&fs::read(dir.join("cut-5"))?)?;
    for ptr in [
        "/frames/5/active/elements/0/directors/0/2",
        "/frames/5/archived/0/elements/0/fields/0/jacobian",
        "/frames/5/active/elements/0/face/offset_m",
        "/frames/5/active/generation",
        "/frames/5/cumulative_applied_work_kn_m",
    ] {
        let mut damaged = saved.clone();
        *damaged.pointer_mut(ptr).unwrap() = serde_json::json!(99);
        fs::write(dir.join("bad"), serde_json::to_vec(&damaged)?)?;
        let failure = cli(&["unexpected", "--resume", "bad"]);
        assert!(!failure.status.success());
        assert!(!dir.join("unexpected").exists());
    }
    fs::remove_dir_all(dir)?;
    Ok(())
}

#[test]
fn five_dof_supported_control_compares_physical_energy_and_response() -> Result<()> {
    use nalgebra::{DMatrix, DVector};
    let geo = patch(0.2, -0.1, 0.);
    let free = [2, 4, 5];
    let mut lift = DMatrix::zeros(36, 15);
    for (j, i) in free.into_iter().enumerate() {
        let d = geo.directors[i];
        let a = V::y().cross(&d).normalize();
        let b = d.cross(&a);
        for axis in 0..3 {
            lift[(i * 6 + axis, j * 5 + axis)] = 1.;
            lift[(i * 6 + 3 + axis, j * 5 + 3)] = a[axis];
            lift[(i * 6 + 3 + axis, j * 5 + 4)] = b[axis];
        }
    }
    let rows: Vec<_> = free.iter().flat_map(|i| 6 * i..6 * i + 6).collect();
    let mut rhs = DVector::zeros(36);
    rhs[6 * 2 + 2] = 1.;
    for ratio in [0.1, 0.01, 0.001] {
        for alpha in [1e-5, 1e-4, 1e-3] {
            let k = Kernel::new(
                &geo,
                &Section {
                    thickness_m: ratio,
                    drilling_factor: alpha,
                    ..section()
                },
                Tying::default(),
            )?;
            let physical = DMatrix::from_column_slice(36, 36, (k.stiffness - k.drill).as_slice());
            let full = DMatrix::from_column_slice(36, 36, k.stiffness.as_slice());
            let reduced = lift.transpose() * &physical * &lift;
            let q5 = &lift
                * reduced
                    .cholesky()
                    .unwrap()
                    .solve(&(lift.transpose() * &rhs));
            let supported = full.select_rows(&rows).select_columns(&rows);
            let q = supported
                .clone()
                .cholesky()
                .unwrap()
                .solve(&rhs.select_rows(&rows));
            let mut q6 = Q36::zeros();
            for (j, &r) in rows.iter().enumerate() {
                q6[r] = q[j];
            }
            let residual = (&supported * &q - rhs.select_rows(&rows)).norm();
            let drill = 0.5 * q6.dot(&(k.drill * q6));
            let p6 = 0.5 * q6.dot(&(k.stiffness * q6)) - drill;
            let p5 = 0.5 * q5.dot(&(&physical * &q5));
            println!(
                "t/L={ratio} alpha={alpha} uz5={} uz6={} physical5={p5} physical6={p6} drilling={drill} residual={residual}",
                q5[14], q6[14]
            );
            assert!(residual < 1e-8 && p5 > 0. && p6 > 0. && drill >= 0.);
        }
    }
    Ok(())
}

#[test]
fn specified_surface_load_preserves_resultant_and_moment() -> Result<()> {
    use soil_fem_3d::curved_shell::surface_action;
    let g = patch(0.2, 0.15, -0.05);
    let top = surface_action(&g, 0.1, 0.05, V::x(), 0.)?;
    let bottom = surface_action(&g, 0.1, -0.05, V::x(), 0.)?;
    let force = |v: Q36| (0..6).map(|i| v[6 * i]).sum::<f64>();
    assert!((force(top) - force(bottom)).abs() > 1e-3);
    // Independent reference from GL product quadrature and finite differences
    // of positions (no rotation operator or analytic director derivatives).
    let points = [
        0.046910077030668,
        0.230765344947158,
        0.5,
        0.769234655052841,
        0.953089922969332,
    ];
    let weights = [
        0.118463442528095,
        0.239314335249683,
        0.284444444444444,
        0.239314335249683,
        0.118463442528095,
    ];
    let omega = V::new(0.21, -0.14, 0.32);
    let q = dofs(&g, |x| (omega.cross(&x), omega));
    let position = |p| -> Result<V> {
        let s = g.surface(p)?;
        Ok(geometry::interp(&s.n, &g.x) + 0.1 * s.d)
    };
    let mut work = 0.;
    for (i, r) in points.into_iter().enumerate() {
        for (j, v) in points.into_iter().enumerate() {
            let s = (1. - r) * v;
            let h = 1e-5;
            let a = (position([r + h, s])? - position([r - h, s])?) / (2. * h);
            let b = (position([r, s + h])? - position([r, s - h])?) / (2. * h);
            work += weights[i]
                * weights[j]
                * (1. - r)
                * a.cross(&b).norm()
                * omega.cross(&position([r, s])?).x;
        }
    }
    assert!((work - top.dot(&q)).abs() < 1e-6);
    assert!(surface_action(&g, 0.1, 0., V::zeros(), 1.).is_err());
    Ok(())
}

#[test]
fn cyclic_numbering_and_world_rotation_preserve_the_curved_element() -> Result<()> {
    let original = patch(0.21, -0.12, 0.02);
    let world = nalgebra::Rotation3::from_euler_angles(0.371, -0.612, 1.193);
    let permutation = [1, 2, 0, 4, 5, 3];
    let rotated = Geometry {
        x: permutation.map(|i| world * original.x[i]),
        directors: permutation.map(|i| world * original.directors[i]),
        e: original.e,
    };
    let a = Kernel::new(&original, &section(), Tying::default())?;
    let b = Kernel::new(&rotated, &section(), Tying::default())?;
    let mut transform = soil_fem_3d::curved_shell::K36::zeros();
    for (j, i) in permutation.into_iter().enumerate() {
        for k in [0, 3] {
            transform
                .fixed_view_mut::<3, 3>(6 * j + k, 6 * i + k)
                .copy_from(world.matrix());
        }
    }
    let diff = transform.transpose() * b.stiffness * transform - a.stiffness;
    assert!(diff.norm() / a.stiffness.norm() < 2e-12);
    Ok(())
}
