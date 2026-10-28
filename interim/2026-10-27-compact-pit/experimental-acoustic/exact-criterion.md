# Frozen continuum MCC branch: all-direction criterion

This additional diagnostic does not determine the actual loading branch and is
not a test of the discrete mixed FEM or its algorithmic tangent. It does not
change any solved state, material parameter or acceptance tolerance.

Let A be the symmetric tensor corresponding to C_elastic : dPhi/dsigma, and D
the positive plastic consistency denominator. For a unit normal n,

Q_elastic = G I + (K+G/3) n n^T,
Q_plastic = Q_elastic - (A n)(A n)^T / D.

The rank-one positive-definiteness criterion is

R(n) = [(n^T A^2 n) - beta (n^T A n)^2] / (G D) < 1,
beta = (K+G/3)/(K+4G/3).

In the principal basis of A, w_i=n_i^2 are nonnegative and sum to one. For a
fixed mean x=sum(a_i w_i), convexity of a_i^2 makes its largest second moment
the secant between a_min and a_max. Hence the global maximum over all normals
reduces to the maximum of

f(x) = (a_min+a_max)x - a_min a_max - beta x^2

on [a_min,a_max]. The maximizing x is the clipped vertex of this concave
quadratic. A normal in the plane of the two extreme principal axes attains it.
The implementation evaluates f as variance+(1-beta)*mean^2 to reduce
cancellation, and independently checks the witness with a full3x3 solve.

For R_max>1, Q_elastic^{-1} A n gives a negative Rayleigh witness. For
R_max<1, every continuum acoustic matrix on this frozen plastic branch is
positive definite. Values within1e-10 of one are reported as numerically
indeterminate. This floating-point classification is not an interval proof.
The witness eigenvalue is not claimed to be the minimum eigenvalue over all
directions. The reduction above is derived here from the rank-one criterion;
it is not attributed to a constitutive-model paper or a commercial solver.
