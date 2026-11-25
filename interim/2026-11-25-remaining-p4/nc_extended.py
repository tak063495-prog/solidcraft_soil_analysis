"""Physical kPa assembly in extended exponent range from homogeneous MCC updates."""
import inspect
from pathlib import Path
import sys
import numpy as np
import nc_restart

Base=nc_restart.Section
source=inspect.getsource(Base.evaluate)
import textwrap
source=textwrap.dedent(source)
# Preserve the entire partition policy, quadrature and original force equations.
needle='        f=-np.einsum'
assert needle.strip() in source
source=source.replace(needle.strip(),'''exp=np.array([s["stress_exponent"] for s in states]).reshape(48,27)
    with np.errstate(under="raise",over="raise",invalid="raise"):
        try:
            stress=np.ldexp(stress.astype(np.longdouble),exp[:,:,None])
            c=np.ldexp(c.astype(np.longdouble),exp[:,:,None,None])
        except FloatingPointError as error:
            raise ArithmeticError("physical extended stress/tangent exponent exhausted") from error
    f=-np.einsum''')
# Accumulate before the single round to the global solver's precision.
source=source.replace('force=np.zeros_like(flat)','force=np.zeros_like(flat,dtype=np.longdouble)').replace(
    'k=np.zeros((len(flat),len(flat)))','k=np.zeros((len(flat),len(flat)),dtype=np.longdouble)').replace(
    'return force[self.free],k[np.ix_(self.free,self.free)]','return force[self.free].astype(float),k[np.ix_(self.free,self.free)].astype(float)')
namespace=dict(Base.evaluate.__globals__)
exec(compile(source,'extended-assembly','exec'),namespace)
class Extended(Base):
    evaluate=namespace['evaluate']

if __name__=='__main__':
    nc_restart.Section=Extended
    options = {} if len(sys.argv)==5 else dict(cap=int(sys.argv[5]),max_ds=float(sys.argv[6]))
    nc_restart.run(Path(sys.argv[1]),Path(sys.argv[2]).resolve(),Path(sys.argv[3]),sys.argv[4],**options)
