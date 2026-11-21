"""Free vibration in SI units. Local Newmark and consistent-mass axial rod FE."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh


def exact(t, omega, zeta=.03, u0=.01, v0=.02):
    if not (omega > 0 and 0 <= zeta < 1):
        raise ValueError("underdamped oscillator required")
    t = np.asarray(t)
    alpha = zeta*omega
    wd = omega*np.sqrt(1-zeta*zeta)
    b = (v0+alpha*u0)/wd
    u = np.exp(-alpha*t)*(u0*np.cos(wd*t)+b*np.sin(wd*t))
    v = np.exp(-alpha*t)*((-alpha*u0+wd*b)*np.cos(wd*t)
                           +(-alpha*b-wd*u0)*np.sin(wd*t))
    a = -2*alpha*v-omega*omega*u
    return u, v, a


def newmark(M, C, K, u0, v0, dt, steps):
    """beta=1/4,gamma=1/2, no forcing; midpoint damping work is explicit."""
    if not (dt > 0 and steps >= 1):
        raise ValueError("positive timestep/steps required")
    n = len(u0)
    u, v, a = [np.zeros((steps+1,n)) for _ in range(3)]
    u[0],v[0] = u0,v0
    a[0] = np.linalg.solve(M,-C@v0-K@u0)
    factor = cho_factor(M + dt*C/2 + dt*dt*K/4)
    dissipation = np.zeros(steps+1)
    for i in range(steps):
        up = u[i]+dt*v[i]+dt*dt*a[i]/4
        vp = v[i]+dt*a[i]/2
        a[i+1] = cho_solve(factor, -C@vp-K@up)
        u[i+1] = up+dt*dt*a[i+1]/4
        v[i+1] = vp+dt*a[i+1]/2
        vm = (v[i+1]+v[i])/2
        dissipation[i+1] = dissipation[i]+dt*(vm@C@vm)
    energy = np.einsum('ni,ij,nj->n',v,M,v)/2 + np.einsum('ni,ij,nj->n',u,K,u)/2
    return {"time_s":np.arange(steps+1)*dt,"u":u,"v":v,"a":a,
            "energy_J":energy,"dissipation_J":dissipation,
            "energy_balance_relative":float(np.max(abs(energy+dissipation-energy[0]))/energy[0])}


def rod(n=8, E=2e9, rho=2000., A=.01, L=2.):
    """Longitudinal rod, u(0)=0, traction(L)=0; no bending model implied."""
    h = L/n
    K, M = np.zeros((n+1,n+1)), np.zeros((n+1,n+1))
    for i in range(n):
        ids = np.ix_([i,i+1],[i,i+1])
        K[ids] += E*A/h*np.array([[1,-1],[-1,1]])
        M[ids] += rho*A*h/6*np.array([[2,1],[1,2]])
    K, M = K[1:,1:], M[1:,1:]
    eigenvalues, vectors = eigh(K,M)
    mode = vectors[:,0]/vectors[-1,0]
    return M,K,mode,float(np.sqrt(eigenvalues[0]))


def diagnose(rod_n=16, periods=5):
    omega = np.sqrt(400/2)
    rows = []
    for zeta in [0.,.03]:
        for per in [20,40,80,160]:
            dt = 2*np.pi/omega/per
            result = newmark(np.array([[2.]]),np.array([[2*zeta*omega*2]]),
                             np.array([[400.]]),np.array([.01]),np.array([.02]),dt,periods*per)
            analytical = exact(result["time_s"],omega,zeta)
            row = {"case":"sdof","zeta":zeta,"steps_per_period":per,"dt_s":dt,
                   "frequency_Hz":omega/(2*np.pi),"undamped_period_s":2*np.pi/omega,
                   "damped_period_s":2*np.pi/(omega*np.sqrt(1-zeta*zeta)),
                   "energy_balance_relative":result["energy_balance_relative"],
                   "initial_energy_J":float(result["energy_J"][0]),
                   "final_energy_J":float(result["energy_J"][-1]),
                   "damping_work_J":float(result["dissipation_J"][-1])}
            for key,truth in zip(["u","v","a"],analytical):
                row[key+"_relative_L2"] = float(np.linalg.norm(result[key][:,0]-truth)/np.linalg.norm(truth))
            rows.append(row)
    spatial = []
    continuous_omega = np.pi/(2*2.)*np.sqrt(2e9/2000.)
    for n in [4,8,16,32]:
        _,_,_,w = rod(n)
        spatial.append({"n":n,"frequency_Hz":w/(2*np.pi),
                        "relative_error":abs(w/continuous_omega-1)})
    M,K,mode,w = rod(rod_n)
    histories = {}
    for per in [20,40,80,160]:
        dt = 2*np.pi/w/per
        r = newmark(M,2*.03*continuous_omega*M,K,.01*mode,.02*mode,dt,periods*per)
        truth = exact(r["time_s"],w,zeta=.03*continuous_omega/w)
        rows.append({"case":"rod_semidiscrete_first_mode","n":rod_n,"steps_per_period":per,
                     "dt_s":dt,"frequency_Hz":w/(2*np.pi),
                     "semidiscrete_damping_ratio":.03*continuous_omega/w,
                     "u_relative_L2":float(np.linalg.norm(r["u"][:,-1]-truth[0])/np.linalg.norm(truth[0])),
                     "v_relative_L2":float(np.linalg.norm(r["v"][:,-1]-truth[1])/np.linalg.norm(truth[1])),
                     "a_relative_L2":float(np.linalg.norm(r["a"][:,-1]-truth[2])/np.linalg.norm(truth[2])),
                     "energy_balance_relative":r["energy_balance_relative"]})
        if per == 160:
            for case,omega_case,rr in [("rod",w,r), ("sdof",omega,newmark(
                    np.array([[2.]]),np.array([[2*.03*omega*2]]),np.array([[400.]]),
                    np.array([.01]),np.array([.02]),2*np.pi/omega/per,periods*per))]:
                histories[case] = {"frequency_Hz":omega_case/(2*np.pi),
                                   "time_s":rr["time_s"].tolist(),
                                   "displacement_m":rr["u"][:,-1].tolist(),
                                   "velocity_m_s":rr["v"][:,-1].tolist(),
                                   "acceleration_m_s2":rr["a"][:,-1].tolist(),
                                   "energy_J":rr["energy_J"].tolist(),
                                   "dissipation_J":rr["dissipation_J"].tolist()}
    for group in [[r for r in rows if r["case"]=="sdof" and r["zeta"]==z] for z in [0.,.03]] + [rows[-4:]]:
        for prev,row in zip(group,group[1:]):
            for key in ["u","v","a"]:
                row["observed_time_order_"+key] = float(np.log2(prev[key+"_relative_L2"]/row[key+"_relative_L2"]))
    for prev,row in zip(spatial,spatial[1:]):
        row["observed_space_order"] = float(np.log2(prev["relative_error"]/row["relative_error"]))
    return {"sdof_and_rod_time":rows,"rod_space":spatial,
            "rod_continuous_frequency_Hz":continuous_omega/(2*np.pi),"histories":histories,
            "accepted":None}
