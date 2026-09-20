"""Serial 1-D adaptive GK15 quadrature matching the legacy cubature rule.

Port of the rule/error estimator and scalar batch-refinement policy in Steven
G. Johnson's cubature hcubature.c (GPL-2.0-or-later; repository LICENSE).
Copyright (c) 2005-2013 Steven G. Johnson. Rule portions based on GNU GSL,
copyright (c) 1996-2000 Brian Gough. Python adaptation: 2026.
Redistribution/modification permitted under GPL version 2 or later.
Provided WITHOUT ANY WARRANTY, including MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the repository LICENSE for the full terms.
The C callback implementation is NOT used. Batched refinement is serial here,
not OpenMP or multiprocessing. Retains the 100000 evaluation ceiling.
"""
import heapq
import math
import sys
import warnings

X = (.9914553711208126392, .9491079123427585245, .8648644233597690728,
     .7415311855993944399, .5860872354676911303, .4058451513773971669,
     .2077849550078984676)
WG = (.1294849661688696933, .2797053914892766679,
      .3818300505051189449, .4179591836734693878)
WK = (.022935322010529225, .0630920926299785533, .1047900103222501838,
      .1406532597155259187, .1690047266392679028, .1903505780647854099,
      .2044329400752988924, .2094821410847278280)
ORDER = (1, 3, 5, 0, 2, 4, 6)

def _rule(f, center, half):
    """Evaluate one embedded GK15 interval and its rescaled absolute error estimate."""
    v0 = f(center)
    pairs = [(f(center-half*X[j]), f(center+half*X[j])) for j in ORDER]
    gauss, kronrod = v0*WG[3], v0*WK[7]
    absolute = abs(kronrod)
    for k, j in enumerate(ORDER):
        a, b = pairs[k]
        if k < 3:
            gauss += WG[k]*(a+b)
        kronrod += WK[j]*(a+b)
        absolute += WK[j]*(abs(a)+abs(b))
    # The absolute-deviation estimate tempers an accidentally small Gauss/Kronrod
    # difference; a machine-precision floor prevents claiming impossible accuracy.
    mean = kronrod*.5
    asc = WK[7]*abs(v0-mean)
    for k, j in enumerate(ORDER):
        a, b = pairs[k]
        asc += WK[j]*(abs(a-mean)+abs(b-mean))
    error = abs(kronrod-gauss)*half
    absolute *= half
    asc *= half
    if asc and error:
        error = asc*min(1., (200*error/asc)**1.5)
    if absolute > sys.float_info.min/(50*sys.float_info.epsilon):
        error = max(error, 50*sys.float_info.epsilon*absolute)
    value = kronrod*half
    if not math.isfinite(value) or not math.isfinite(error):
        raise RuntimeError("Non-finite quadrature result")
    return value, error

def integrate(f, low, high, rtol=1e-8, maxeval=100000):
    """Adaptively integrate a scalar callable over finite increasing bounds.

    Refine the largest-error intervals in batches, but execute serially.
    At maxeval return the current estimate with a warning, as documented;
    this is not evidence that the requested relative tolerance was reached.
    """
    if not math.isfinite(rtol) or rtol <= 0 or maxeval < 15:
        raise ValueError("Require positive finite tolerance and at least 15 evaluations")
    if not math.isfinite(low+high) or low >= high:
        raise ValueError("Quadrature bounds must be finite and increasing")
    center, half = .5*(high+low), .5*(high-low)
    total, error = _rule(f, center, half)
    serial = 0
    heap = [(-error, serial, center, half, total)]
    evaluations = 15
    while evaluations < maxeval and error > abs(total)*rtol:
        remaining_error = error
        original_total = total
        regions = []
        while heap:
            negative_err, _, center, half, value = heapq.heappop(heap)
            e = -negative_err
            remaining_error -= e
            total -= value
            error -= e
            half *= .5
            regions.extend(((center-half, half), (center+half, half)))
            evaluations += 30
            if remaining_error <= abs(original_total)*rtol or not heap or evaluations >= maxeval:
                break
        for center, half in regions:
            value, e = _rule(f, center, half)
            total += value
            error += e
            serial += 1
            heapq.heappush(heap, (-e, serial, center, half, value))
    if error > abs(total)*rtol:
        warnings.warn("Quadrature evaluation limit reached; legacy estimate returned",
                      RuntimeWarning, stacklevel=2)
    return total
