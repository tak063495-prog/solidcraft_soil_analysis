"""Unit-area, unit-height shear column in series with a friction interface.

Three tangential nodes: ground(0), interface/column(1), prescribed top(2).
Normal column strain is prescribed; a prescribed compressive clamp acts at the
interface, with zero normal gap/work. Fixed hydrostatic pore pressure. This is
a diagnostic reduced model, not the original 12-Tet problem.
"""
import numpy as np
from scipy.optimize import brentq
from mc import Material, implicit, hydrostatic


def solve(material_parts=4, contact_parts=4, steps=24):
    if min(material_parts, contact_parts, steps) < 1:
        raise ValueError('positive partitions required')
    mat = Material()
    kt, mu, p = 200000., .3, hydrostatic(0.)
    # Start stress-free, then apply normal strain and shear through full history.
    s, ep = np.zeros((3, 3)), np.zeros((3, 3))
    u1, u2, epsn, q, slip, clamp = 0., 0., 0., 0., 0., 200.
    U0 = 0.
    W, Db, Df = 0., 0., 0.
    history = []

    def prescribed(t):
        return (.012*np.sin(2*np.pi*t),
                -.0025*(1-np.cos(np.pi*t))-.0006*np.sin(4*np.pi*t),
                200.+150.*np.sin(2*np.pi*t+.4))

    # Clamp at t=0 must match initial state; no tangential work at startup.
    clamp = prescribed(0.)[2]
    for k in range(1, steps+1):
        t0, t1 = (k-1)/steps, k/steps
        top, en, new_clamp = prescribed(t1)

        def evaluate(v):
            bulk = s.copy()
            dep = np.zeros((3, 3))
            wb, wn, db, min_db = 0., 0., 0., 0.
            for j in range(material_parts):
                a, b = j/material_parts, (j+1)/material_parts
                # Resolved physical normal strain path; unknown middle-node
                # tangential displacement is interpolated between endpoints.
                ea = prescribed(t0+(t1-t0)*a)[1]
                eb = prescribed(t0+(t1-t0)*b)[1]
                de = np.diag([0., 0., eb-ea])
                de[0, 2] = de[2, 0] = ((top-u2)-(v-u1))/(2*material_parts)
                r = implicit(bulk, de, mat)
                midpoint = .5*(bulk+r['stress'])
                d = float(np.sum(midpoint*r['dep']))
                wb += float(np.sum(midpoint*de))
                wn += float(midpoint[2,2]*de[2,2])
                db += d
                min_db = min(min_db, d)
                dep += r['dep']
                bulk = r['stress']
            traction, ds, wf, df, min_df = q, 0., 0., 0., 0.
            for j in range(contact_parts):
                b = (j+1)/contact_parts
                effective_n = max(prescribed(t0+(t1-t0)*b)[2]-p, 0.)
                du = (v-u1)/contact_parts
                trial = traction+kt*du
                new_q = float(np.clip(trial, -mu*effective_n, mu*effective_n))
                dslip = du-(new_q-traction)/kt
                # Midpoint conjugate work gives exact spring energy identity.
                diss = .5*(traction+new_q)*dslip
                wf += .5*(traction+new_q)*du
                df += diss
                min_df = min(min_df, diss)
                ds += dslip
                traction = new_q
            return dict(stress=bulk, dep=dep, traction=traction, dslip=ds,
                        wb=wb, wn=wn, db=db, wf=wf, df=df, min_db=min_db, min_df=min_df)

        # Solve all free-node tangential equilibrium at each accepted endpoint.
        def residual(v):
            r = evaluate(v)
            return r['stress'][0, 2]-r['traction']
        v = brentq(residual, -1., 1., xtol=1e-14)
        r = evaluate(v)
        # Independent external reaction quadrature on the fixed global grid.
        W += .5*(s[0, 2]+r['stress'][0, 2])*(top-u2)
        # Total stress includes -p I: fixed pressure work is tracked separately.
        W += r['wn']
        water_work = -p*(en-epsn)
        old_U = .5*float(np.sum(s*mat.compliance(s)))+.5*q*q/kt
        s, ep = r['stress'], ep+r['dep']
        q, slip = r['traction'], slip+r['dslip']
        Db, Df = Db+r['db'], Df+r['df']
        U = .5*float(np.sum(s*mat.compliance(s)))+.5*q*q/kt
        # Forces [ground reaction, internal residual, top reaction] with all
        # reactions included; no omitted constrained endpoint in balance check.
        nodal = np.array([-q, q-s[0, 2], s[0, 2]])
        total_zz = s[2,2]-p
        vectors = np.array([[-q,0.,-new_clamp],
                            [q-s[0,2],0.,new_clamp-total_zz],
                            [s[0,2],0.,total_zz]])
        history.append(dict(stage='cycle', increment=k, time=t1,
                            u_middle=v, u_top=top, stress=s.tolist(),
                            plastic_strain=ep.tolist(), plastic_slip=slip,
                            traction=q, total_clamp=new_clamp,
                            effective_contact_pressure=max(new_clamp-p,0.),
                            pore=p, nodal_forces=nodal.tolist(),
                            all_nodal_force_vectors=vectors.tolist(),
                            free_force_residual=abs(nodal[1]),
                            all_node_force_sum=float(np.linalg.norm(vectors.sum(axis=0))),
                            elastic_energy=U, plastic_dissipation=Db,
                            friction_dissipation=Df, external_skeleton_work=W,
                            cumulative_pressure_work=-p*en,
                            external_total_work=W-p*en,
                            energy_defect=W-(U-U0)-Db-Df,
                            local_energy_identity_error=r['wb']+r['wf']-(U-old_U)-r['db']-r['df'],
                            minimum_bulk_dissipation_increment=r['min_db'],
                            minimum_friction_dissipation_increment=r['min_df'],
                            local_internal_work=r['wb']+r['wf'],
                            pressure_work_increment=water_work))
        u1, u2, epsn, clamp = v, top, en, new_clamp
    return history
