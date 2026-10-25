use nalgebra::{Matrix3, Vector3};
use soil_fem_3d::{
    Result,
    embedding::{Binding, Point},
    mesh::Mesh,
};

fn mesh() -> Result<Mesh> {
    Ok(serde_json::from_str(include_str!(
        "../examples/mc-column-mesh.json"
    ))?)
}

#[test]
fn tet4_and_tet10_embedding_preserve_affine_fields_force_moment_and_work() -> Result<()> {
    for quadratic in [false, true] {
        let linear = mesh()?;
        let mesh = if quadratic {
            linear.elevate_straight()?
        } else {
            linear
        };
        let p = Point {
            position_m: [0.23, 0.31, 0.67],
            body: None,
        };
        let b = Binding::locate(&mesh, &p, 1e-10)?;
        let a = Matrix3::new(0.1, -0.2, 0.3, 0.15, 0.13, -0.17, -0.07, 0.04, 0.12);
        let offset = Vector3::new(0.03, -0.02, 0.01);
        let u: Vec<[f64; 3]> = mesh
            .nodes
            .iter()
            .map(|x| (a * Vector3::from(*x) + offset).into())
            .collect();
        assert!((b.displacement(&u) - (a * Vector3::from(p.position_m) + offset)).norm() < 1e-11);
        if quadratic {
            assert!(b.weights.iter().any(|(_, w)| *w < 0.));
        }
        let f = Vector3::new(2., -3., 4.);
        let mut nodal = vec![[0.; 3]; mesh.nodes.len()];
        b.distribute(f, &mut nodal);
        let sum = nodal
            .iter()
            .fold(Vector3::zeros(), |s, v| s + Vector3::from(*v));
        let moment = mesh
            .nodes
            .iter()
            .zip(&nodal)
            .fold(Vector3::zeros(), |s, (x, v)| {
                s + Vector3::from(*x).cross(&Vector3::from(*v))
            });
        let work: f64 = u
            .iter()
            .zip(&nodal)
            .map(|(u, f)| Vector3::from(*u).dot(&Vector3::from(*f)))
            .sum();
        assert!((sum - f).norm() < 1e-12);
        assert!((moment - Vector3::from(p.position_m).cross(&f)).norm() < 1e-10);
        assert!((work - b.displacement(&u).dot(&f)).abs() < 1e-12);
    }
    Ok(())
}

#[test]
fn curved_tet10_inverse_and_interface_hosts_keep_coordinate_and_activity() -> Result<()> {
    let mut mesh = mesh()?.elevate_straight()?;
    for p in &mut mesh.nodes {
        p[0] += 0.15 * p[1] * (1. - p[1]);
    }
    mesh.validate()?;
    for mut x in [[0.23, 0.31, 0.67], [0.5, 0.5, 0.5], [0., 0., 0.]] {
        x[0] += 0.15 * x[1] * (1. - x[1]);
        let b = Binding::locate(
            &mesh,
            &Point {
                position_m: x,
                body: None,
            },
            1e-11,
        )?;
        assert!(b.maximum_coordinate_residual_m < 1e-11);
        let mut active = vec![false; mesh.tetrahedra.len()];
        assert!(!b.is_active(&active));
        active[*b.host_elements.last().unwrap()] = true;
        assert!(b.is_active(&active));
    }
    Ok(())
}

#[test]
fn embedding_rejects_missing_host_body_outside_and_nonconforming_overlap() -> Result<()> {
    let mut mesh = mesh()?;
    assert!(
        Binding::locate(
            &mesh,
            &Point {
                position_m: [0.2; 3],
                body: Some("missing".into())
            },
            1e-8
        )
        .is_err()
    );
    assert!(
        Binding::locate(
            &mesh,
            &Point {
                position_m: [3.; 3],
                body: None
            },
            1e-8
        )
        .is_err()
    );
    let tet = mesh.tetrahedra[0];
    let center = std::array::from_fn(|a| tet.iter().map(|&n| mesh.nodes[n][a]).sum::<f64>() / 4.);
    let start = mesh.nodes.len();
    let duplicates: Vec<_> = tet.iter().map(|&n| mesh.nodes[n]).collect();
    mesh.nodes.extend(duplicates);
    mesh.tetrahedra
        .push([start, start + 1, start + 2, start + 3]);
    mesh.material_ids.push(mesh.material_ids[0]);
    mesh.region_ids.push(mesh.region_ids[0]);
    assert!(
        Binding::locate(
            &mesh,
            &Point {
                position_m: center,
                body: None
            },
            1e-10
        )
        .is_err()
    );
    Ok(())
}
