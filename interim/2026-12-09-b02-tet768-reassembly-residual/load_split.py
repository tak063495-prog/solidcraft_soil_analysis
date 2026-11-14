"""Discarded linear RHS source probes on the same independently assembled system."""
import numpy as np
from scipy.linalg import lu_solve
import fem_reference as ref

def diagnose_sources(model,state,binding,qb,factorization,mesh,case,active,free):
    # Boundary load derivative only: reference positions and constant gravity
    # stay fixed. The source responses never replace the original history.
    full=model.assemble(state,1.)[2]-model.assemble(state,0.)[2]
    nodal=np.zeros(model.size)
    for p in case['beam_system']['nodal_loads']:
        n=model.node_ids[p['node_id']]
        nodal[6*n:6*n+6]+=np.r_[p.get('boundary_force_kn',[0.]*3),p.get('boundary_moment_kn_m',[0.]*3)]
    soil=ref.conditions_load(mesh,case['boundary'],active).ravel()[free]
    tip_rhs=np.r_[binding.T@nodal,qb.T@nodal]
    other_rhs=np.r_[soil+binding.T@(full-nodal),qb.T@(full-nodal)]
    response=lu_solve(factorization,np.column_stack((tip_rhs,other_rhs)))
    ds=binding@response[:len(free)]+qb@response[len(free):]
    tip=model.node_ids['tip'];attachment=model.data[0][0][0]
    return dict(kind='same input-only matrix, RHS derivative probes; no adopted state',
                nodal_rhs_compliance_kn_m=float(tip_rhs@response[:,0]),other_rhs_compliance_kn_m=float(other_rhs@response[:,1]),
                cross_compliance_kn_m=float(tip_rhs@response[:,1]),
                nodal_structural_response=ds[:,0].reshape(-1,6).tolist(),other_structural_response=ds[:,1].reshape(-1,6).tolist(),
                nodal_tip_response=ds[6*tip:6*tip+6,0].tolist(),nodal_attachment_response=ds[6*attachment:6*attachment+6,0].tolist(),
                nodal_load_work_kn_m=float(nodal@ds[:,0]),nodal_full_work_kn_m=float(full@ds[:,0]),
                source_probe_linear_solves=2,independent_history_updated=False)
