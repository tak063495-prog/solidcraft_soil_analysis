"""Exact-reference self checks. Passing these does not accept any missing feature.

Source TODO = interim/2026-11-16-issue-triage/todo.ja.md.
(a) = interim/2026-11-13-joint-ext-combos/README.ja.md.
(b) = interim/2026-11-14-plate6-curved-offset/README.ja.md.
Each item cites the TODO and the contract in its test docstring. Numerical
gates below are the existing 1e-10 work/constraint/reassembly gates, not new
acceptance tolerances for shell benchmarks or unsaturated convergence.
"""
import numpy as np
import pytest
from oracles import tri6, triangle_rule, hydrostatic_triangle, darcy, diffusion, strip


@pytest.mark.parametrize('level', [1.5, 2.])
def test_A01_hydrostatic_all_node_force_and_work_reference(level):
    """todo A-01; (a) README: p=gamma(level-z), all-node force/work 1e-10."""
    p = np.array([[0.,0.,0.], [2.,0.,0.], [0.,0.,1.], [1.,0.,0.], [1.,0.,.5], [0.,0.,.5]])
    gamma = 9.81
    f = hydrostatic_triangle(p, level, gamma)
    # Area=1, oriented normal=-y; exact integral of a linear pressure.
    np.testing.assert_allclose(f.sum(axis=0), [0., -gamma*(level-1/3), 0.], atol=1e-10, rtol=0.)
    # Virtual translation v_y=z tests a second pressure moment independently.
    # Integral z=1/3, z^2=1/6 on this triangle.
    assert abs(np.sum(f[:,1]*p[:,2])+gamma*(level/3-1/6)) < 1e-10
    # For constant pressure the corner shape-function integrals vanish.
    horizontal = np.array([[0.,0.,0.],[2.,0.,0.],[0.,1.,0.],[1.,0.,0.],[1.,.5,0.],[0.,.5,0.]])
    constant = hydrostatic_triangle(horizontal, level, gamma)
    np.testing.assert_allclose(constant[:3],0.,atol=1e-10,rtol=0.)
    np.testing.assert_allclose(constant[3:,2],gamma*level/3,atol=1e-10,rtol=0.)


@pytest.mark.parametrize('top,expected_flux', [(2.,0.), (2.5,-5e-7), (3.,-1e-6)])
def test_A03_darcy_head_flux_and_birth_subdomain_reference(top,expected_flux):
    """todo A-03; (a) README: independent water balance; head gate 1e-8 m."""
    z = np.linspace(0., 1., 11)
    h, pressure, q = darcy(z, 2., top, 1., 1e-6, 9.81)
    np.testing.assert_allclose(h[[0,-1]], [2.,top], atol=1e-8, rtol=0.)
    np.testing.assert_allclose(pressure, 9.81*(h-z), atol=1e-10, rtol=0.)
    # Every new subdomain shares the same constant flux, no invented source.
    assert abs(q-expected_flux) < 1e-16  # existing Richards absolute-flow gate
    for lower, upper in [(0.,.5), (.5,1.)]:
        hh, _, qq = darcy(np.array([lower,upper]), 2.,top,1.,1e-6,9.81)
        assert abs(qq+1e-6*(hh[1]-hh[0])/(upper-lower)) < 1e-10


@pytest.mark.parametrize('time', [0., 100., 1000.])
def test_A04_saturated_richards_transient_storage_balance_reference(time):
    """todo A-04; saturated manufactured subcase; NOT an unsaturated JOINT pass.

    (a) README requires transient storage and independent mass balance. This
    checks the identity analytically, without inventing a convergence gate.
    """
    length, k, storage, amplitude = 1., 1e-6, 1e-3, .1
    z, weights = np.polynomial.legendre.leggauss(20)
    z, weights = (z+1)/2, weights/2
    h, ht, _ = diffusion(z,time,length,k,storage,amplitude)
    _, _, ends = diffusion(np.array([0.,length]),time,length,k,storage,amplitude)
    assert abs(storage*(weights@ht)-(ends[0]-ends[1])) < 1e-10
    assert np.all(2.+h-z > 0.)
    boundary, _, _ = diffusion(np.array([0.,length]),time,length,k,storage,amplitude)
    np.testing.assert_allclose(boundary,0.,atol=1e-8,rtol=0.)


@pytest.mark.parametrize('kind', ['cylinder','sphere'])
def test_B01_curved_tri6_metric_and_rigid_motion_reference(kind):
    """todo B-01; (b) README x=sum(N*x), a_r/a_s, geometric normals.

    This is a quadratic interpolant of shell nodes, not exact circle geometry
    and not a Scordelis-Lo/pinch/hemisphere benchmark result.
    """
    rs = np.array([[0.,0.],[.2,0.],[0.,.2],[.1,0.],[.1,.1],[0.,.1]])
    r,s = rs.T
    points = (np.c_[np.sin(r),s,np.cos(r)] if kind=='cylinder' else
              np.c_[np.sin(r)*np.cos(s),np.sin(s),np.cos(r)*np.cos(s)])
    omega = np.array([.2,-.3,.1]); translation = np.array([.1,.2,.3])
    displacements = translation+np.cross(omega,points)
    normals=[]
    for r,s,_ in triangle_rule():
        n,dr,ds = tri6(r,s)
        ar, ass = dr@points, ds@points
        dur, dus = dr@displacements, ds@displacements
        metric_change = [2*ar@dur, ar@dus+ass@dur, 2*ass@dus]
        np.testing.assert_allclose(metric_change,0.,atol=1e-10,rtol=0.)
        normal = np.cross(ar,ass); assert np.linalg.norm(normal)>0.
        normal /= np.linalg.norm(normal); normals.append(normal)
        np.testing.assert_allclose(n@displacements, translation+np.cross(omega,n@points),atol=1e-10,rtol=0.)
    assert np.linalg.norm(normals[0]-normals[-1])>1e-3


@pytest.mark.parametrize('support', ['cantilever','simple'])
@pytest.mark.parametrize('offset', [-.05,0.,.05])
def test_B02_offset_strip_exact_reference(support, offset):
    """todo B-02; (b) README analytic strip, original structural 1e-9 gate."""
    young, thickness, force, length = 2e8,.1,10.,2.
    x=np.array([0.,1.,length])
    u,w,theta,curvature=strip(x,young,thickness,offset,force,length,support)
    # Fixed small reference values in m/rad, not values generated by the
    # function under test: E=2e8 kPa,t=.1 m,N=10 kN/m,L=2 m.
    sign=np.sign(offset)
    np.testing.assert_allclose(u,([0.,5e-7,1e-6] if offset==0. else [0.,2e-6,4e-6]),atol=1e-9,rtol=0.)
    if support=='cantilever':
        np.testing.assert_allclose(w,sign*np.array([0.,-1.5e-5,-6e-5]),atol=1e-9,rtol=0.)
        np.testing.assert_allclose(theta,sign*np.array([0.,3e-5,6e-5]),atol=1e-9,rtol=0.)
    else:
        np.testing.assert_allclose(w,sign*np.array([0.,1.5e-5,0.]),atol=1e-9,rtol=0.)
        np.testing.assert_allclose(theta,sign*np.array([-3e-5,0.,3e-5]),atol=1e-9,rtol=0.)
    assert abs(young*thickness*(u[-1]/length-offset*curvature)-force)<1e-8
    assert abs(young*thickness**3/12*curvature-offset*force)<1e-8
    if support=='cantilever':
        np.testing.assert_allclose([u[0],w[0],theta[0]],0.,atol=1e-9,rtol=0.)
    else:
        np.testing.assert_allclose(w[[0,-1]],0.,atol=1e-9,rtol=0.)
    if offset==0.:
        np.testing.assert_array_equal(w,np.zeros(3))
