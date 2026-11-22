"""T8 §4 estimates for a signed scalar measured at the same physical location."""
import math


def richardson(values, widths, *, similar, complete, reference_converged,
               other_axes_fixed, smooth):
    """Coarse to fine. The caller identifies the case, observable and physical widths."""
    result = dict(observed_order=None, extrapolated=None, gci_fine=None, holds=[])
    holds = result["holds"]
    for valid, reason in ((similar, "non-similar grids"), (complete, "incomplete series"),
                          (reference_converged, "unresolved reference"),
                          (other_axes_fixed, "mixed refinement axes"), (smooth, "non-smooth observable")):
        if not valid:
            holds.append(reason)
    if len(values) != 3 or len(widths) != 3:
        holds.append("three levels required")
        return result
    if any(q is None or not math.isfinite(q) for q in values):
        holds.append("missing finite observable")
        return result
    if any(not math.isfinite(h) or h <= 0 for h in widths):
        holds.append("invalid physical widths")
        return result
    coarse, middle, fine = values
    r, s = widths[0] / widths[1], widths[1] / widths[2]
    if min(r, s) <= 1 or abs(r - s) > max(math.ulp(r), math.ulp(s)):
        holds.append("unequal refinement ratios")
    a, b = coarse - middle, middle - fine
    if a * b < 0:
        holds.append("signed differences reverse")
    if abs(a) <= max(math.ulp(coarse), math.ulp(middle)) or abs(b) <= max(math.ulp(middle), math.ulp(fine)):
        holds.append("roundoff-scale differences")
    if holds:
        return result
    order = math.log(abs(a / b)) / math.log(r)
    result["observed_order"] = order
    if order <= 0:
        holds.append("non-positive order")
        return result
    correction = (fine - middle) / math.expm1(order * math.log(r))
    result["extrapolated"] = fine + correction
    result["dimensional_correction"] = abs(correction)
    if abs(fine) > max(math.ulp(coarse), math.ulp(middle), math.ulp(fine)):
        result["gci_fine"] = 1.25 * abs(correction / fine)
    return result
