"""Script 04's coherent two-disk model and scalar-diagonal data objective.

All coordinates here are tangent-plane arcsec (east, north); PA is measured
east of north. Node values are nonnegative intrinsic speeds, not LOS speeds.
The two permanent component labels have opposite, fixed rotation signs.
No pPXF fits or noise recalibrations happen in this module.
"""
from __future__ import annotations

from dataclasses import dataclass
import time
import logging

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.signal import fftconvolve


def disk_coordinates(east, north, pa_deg, inclination_deg):
    """Return disk radius and cos(azimuth), with astronomical PA conventions."""
    angle = np.deg2rad(pa_deg)
    major = east * np.sin(angle) + north * np.cos(angle)
    minor = east * np.cos(angle) - north * np.sin(angle)
    radius = np.hypot(major, minor / np.cos(np.deg2rad(inclination_deg)))
    cosine = np.divide(major, radius, out=np.zeros_like(radius), where=radius > 0)
    return radius, cosine


def radial_nodes(first, spacing, target, *, outward):
    """Freeze a node grid before optimization; full coverage rounds outward."""
    if min(first, spacing) <= 0 or not np.all(np.isfinite([first, spacing, target])):
        raise ValueError("Positive finite node scales and a finite extent are required")
    offset = max(0., (target - first) / spacing)
    count = int(np.ceil(offset) if outward else np.floor(offset + 0.5))
    return first + spacing * np.arange(count + 1)


def curvature_penalty(speed_a, speed_b, coefficient):
    """Published squared second differences, excluding the imposed origin.

    The first stencil is (U1,U2,U3), not (0,U1,U2). C=lambda*A has units
    (km/s)^-2 when node speeds are supplied in km/s. Zero means the control.
    """
    if coefficient < 0 or not np.isfinite(coefficient):
        raise ValueError("The effective regularity coefficient must be finite and nonnegative")
    roughness = sum(float(np.sum(np.diff(u, n=2)**2)) for u in (speed_a, speed_b))
    return coefficient * roughness, roughness


def profile_cube(chi2, va_axis, vb_axis, fractions, va, vb):
    """Interpolate EACH fraction plane, then profile; never extrapolate.

    chi2 has shape (bin, VA, VB, fraction). At exact grid coordinates, failed
    corners with zero interpolation weight do not invalidate a valid fit.
    A failed fraction plane is unavailable; at least one plane must be valid.
    """
    va, vb = np.asarray(va), np.asarray(vb)
    nbin = len(va)
    if chi2.shape[0] != nbin:
        raise ValueError("Predictions and likelihood cubes have different bin counts")
    supported = (np.isfinite(va) & np.isfinite(vb) &
                 (va >= va_axis[0]) & (va <= va_axis[-1]) &
                 (vb >= vb_axis[0]) & (vb <= vb_axis[-1]))
    if not np.all(supported):
        return None
    ia = np.clip(np.searchsorted(va_axis, va, side="right") - 1, 0, len(va_axis)-2)
    ib = np.clip(np.searchsorted(vb_axis, vb, side="right") - 1, 0, len(vb_axis)-2)
    ta = (va-va_axis[ia])/(va_axis[ia+1]-va_axis[ia])
    tb = (vb-vb_axis[ib])/(vb_axis[ib+1]-vb_axis[ib])
    scores = np.zeros((nbin, len(fractions)))
    available = np.ones_like(scores, dtype=bool)
    bids = np.arange(nbin)
    for da, db, weights in ((0,0,(1-ta)*(1-tb)), (1,0,ta*(1-tb)),
                            (0,1,(1-ta)*tb), (1,1,ta*tb)):
        corner = chi2[bids, ia+da, ib+db, :]
        finite = np.isfinite(corner)
        active = weights[:, None] > 0
        available &= finite | ~active
        scores += np.where(finite, corner, 0.) * weights[:, None]
    scores[~available] = np.inf
    chosen = np.argmin(scores, axis=1)
    best = scores[bids, chosen]
    if not np.all(np.isfinite(best)):
        return None
    return best, np.asarray(fractions)[chosen]


@dataclass
class DiskProblem:
    """One fixed dataset/extent/PSF configuration, safe to send to workers.

    theta = [east0,north0,PA,inclination,Vsys, U_A[1:K], U_B[1:K]].
    bin_index addresses the full spatial arrays at the actual RH3 members.
    Selection is fixed: unsupported trials fail rather than dropping bins.
    """
    east: np.ndarray
    north: np.ndarray
    light: np.ndarray
    member_y: np.ndarray
    member_x: np.ndarray
    member_bin: np.ndarray
    selected: np.ndarray
    bin_ids: np.ndarray
    chi2: np.ndarray
    va_axis: np.ndarray
    vb_axis: np.ndarray
    fractions: np.ndarray
    nodes: np.ndarray
    bounds: np.ndarray
    initial_geometry: np.ndarray
    initial_speed: float
    sign_a: int
    kernel: np.ndarray | None = None

    @property
    def nnode(self):
        return len(self.nodes)

    def predict(self, theta, *, diagnostics=True):
        """Return bin LOS means/shear and radial-padding diagnostic fractions."""
        x0,y0,pa,inc,vsys = theta[:5]
        radius, cosine = disk_coordinates(self.east-x0, self.north-y0, pa, inc)
        projection = np.sin(np.deg2rad(inc))*cosine
        denominator = self.light
        if self.kernel is not None:
            denominator = fftconvolve(self.light, self.kernel, mode="same")
        y,x,b = self.member_y,self.member_x,self.member_bin
        nb = len(self.bin_ids)
        total = np.bincount(b, weights=denominator[y,x], minlength=nb)
        if np.any(total[self.selected] <= 0):
            return None
        velocities, shears = [], []
        for sign, speeds in ((self.sign_a,theta[5:5+self.nnode]),
                             (-self.sign_a,theta[5+self.nnode:])):
            # np.interp explicitly holds the outermost speed constant. The
            # geometry and PSF padding never create a new fitted radial node.
            speed = np.interp(radius, np.r_[0.,self.nodes], np.r_[0.,speeds])
            field = vsys + sign*speed*projection
            first = self.light*field
            if self.kernel is not None:
                first = fftconvolve(first,self.kernel,mode="same")
            mean = np.divide(np.bincount(b,weights=first[y,x],minlength=nb),total,
                             out=np.full(nb,np.nan),where=total>0)
            velocities.append(mean[self.selected])
            if diagnostics:
                second = self.light*field**2
                if self.kernel is not None:
                    second = fftconvolve(second,self.kernel,mode="same")
                square = np.divide(np.bincount(b,weights=second[y,x],minlength=nb),total,
                                   out=np.full(nb,np.nan),where=total>0)
                shears.append(np.sqrt(np.maximum(0.,square[self.selected]-mean[self.selected]**2)))
        if not diagnostics:
            return (*velocities, None, None, None)
        # Report the fraction of received model light originating beyond the
        # last node, including photons carried INTO an aperture by the PSF.
        padded_light = self.light*(radius > self.nodes[-1])
        if self.kernel is not None:
            padded_light = fftconvolve(padded_light,self.kernel,mode="same")
        padded = np.bincount(b,weights=padded_light[y,x],minlength=nb)
        fraction = np.clip(np.divide(padded,total,out=np.zeros(nb),where=total>0),0.,1.)
        return (*velocities,*shears,fraction[self.selected])

    def evaluate(self, theta, coefficient, *, details=False):
        """Return Q_data+P, or a rejection with no partial-bin objective."""
        theta = np.asarray(theta)
        if (theta.shape != (len(self.bounds),) or not np.all(np.isfinite(theta)) or
                np.any(theta < self.bounds[:,0]) or np.any(theta > self.bounds[:,1])):
            return None if details else np.inf
        prediction = self.predict(theta, diagnostics=details)
        if prediction is None:
            return None if details else np.inf
        fitted = profile_cube(self.chi2,self.va_axis,self.vb_axis,self.fractions,
                              prediction[0],prediction[1])
        if fitted is None:
            return None if details else np.inf
        contributions, fraction = fitted
        penalty, roughness = curvature_penalty(theta[5:5+self.nnode],theta[5+self.nnode:],coefficient)
        score = float(contributions.sum())
        if not details:
            return score+penalty
        return dict(data_score=score,penalty=penalty,roughness=roughness,
                    objective=score+penalty,bin_score=contributions,fraction=fraction,
                    va=prediction[0],vb=prediction[1],shear_a=prediction[2],
                    shear_b=prediction[3],radial_padding_fraction=prediction[4])


def starting_models(problem, *, seed, population_factor):
    """MaNGA-informed shapes, amplitude/asymmetry variations, broad starts.

    These are search initializers, not extra likelihood terms or priors. Keep
    identical seeds/settings across regularity strengths for fair comparisons.
    """
    rng = np.random.default_rng(seed)
    k = problem.nnode
    radius = problem.nodes/problem.nodes[-1]
    shapes = (np.ones(k), np.minimum(1.,radius*3.), 1.-0.3*radius)
    base = problem.initial_speed
    candidates = []
    for shape in shapes:
        for amplitude in (0.5,0.75,1.,1.25,1.5):
            for asymmetry in (0.8,1.,1.25):
                u = base*amplitude*shape
                state = np.r_[problem.initial_geometry,u,u*asymmetry]
                state = np.clip(state,problem.bounds[:,0],problem.bounds[:,1])
                if np.isfinite(problem.evaluate(state,0.)):
                    candidates.append(state)
    if not candidates:
        raise ValueError("No supported initial disk model. Inspect geometry, velocity frame and cube support.")
    # Even very short curves need starts spanning the allowed geometry; do
    # not fill the entire population with identical-geometry MaNGA seeds.
    count = max(len(candidates)+10,population_factor*len(problem.bounds),5)
    attempts = 0
    while len(candidates) < count and attempts < 30*count:
        attempts += 1
        state = candidates[rng.integers(len(candidates))].copy()
        if attempts % 3 == 0:
            state[:5] = rng.uniform(problem.bounds[:5,0],problem.bounds[:5,1])
            state[5:] = rng.uniform(0.,problem.bounds[5:,1])
        else:
            scale = (problem.bounds[:,1]-problem.bounds[:,0])
            state += rng.normal(0.,0.12,len(state))*scale
            state = np.clip(state,problem.bounds[:,0],problem.bounds[:,1])
        if np.isfinite(problem.evaluate(state,0.)):
            candidates.append(state)
    return np.asarray(candidates[:max(count,5)])


def fit_disk(problem, coefficient, *, seed=1, population_factor=6,
             global_iterations=80, local_starts=4, local_iterations=400, tolerance=1e-5):
    """Population exploration, then bounded derivative-free local refinement.

    Differential evolution samples SEARCH states; it does not produce posterior
    draws. Powell refines several promising starts against the SAME objective.
    Non-convergence is saved as a flag, never presented as a proved optimum.
    """
    start = time.monotonic()
    if min(population_factor,local_starts,local_iterations) < 1 or global_iterations < 0:
        raise ValueError('Positive population/local budgets and nonnegative global iterations are required')
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Optimizer tolerance must be finite and positive')
    logger = logging.getLogger(__name__)
    history = []
    rejected = 0
    def objective(theta):
        nonlocal rejected
        value = problem.evaluate(theta,coefficient)
        if not np.isfinite(value):
            rejected += 1
        return value
    population = starting_models(problem,seed=seed,population_factor=population_factor)
    def callback(theta, convergence):
        result = problem.evaluate(theta,coefficient,details=True)
        if result is None:
            return False
        history.append(dict(seconds=time.monotonic()-start,objective=result['objective'],
                            data_score=result['data_score'],penalty=result['penalty']))
        if len(history)==1 or len(history)%5==0:
            logger.info('step=%d Qdata=%.8g P=%.8g geometry=%s elapsed=%.1fs',
                        len(history),result['data_score'],result['penalty'],
                        np.asarray(theta[:5]).round(5).tolist(),time.monotonic()-start)
        return False
    explored = differential_evolution(objective,problem.bounds,init=population,
                                     maxiter=global_iterations,tol=tolerance,seed=seed,
                                     polish=False,callback=callback,workers=1)
    pool = np.vstack((explored.x,explored.population,population))
    scores = np.array([objective(p) for p in pool])
    order = np.argsort(scores)
    starts = []
    span = np.maximum(problem.bounds[:,1]-problem.bounds[:,0],1.)
    for index in order:
        if np.isfinite(scores[index]) and all(np.linalg.norm((pool[index]-p)/span)>1e-3 for p in starts):
            starts.append(pool[index])
        if len(starts) >= local_starts:
            break
    alternatives = []
    for state in starts:
        refined = minimize(objective,state,method='Powell',bounds=problem.bounds,
                           callback=lambda theta:callback(theta,None),
                           options=dict(maxiter=local_iterations,ftol=tolerance,xtol=tolerance))
        # Powell can terminate at an unsupported point on a constrained surface.
        # Preserve the valid starting state instead of accepting a worse result.
        accepted = refined.x if np.isfinite(refined.fun) and refined.fun <= objective(state) else state
        detail = problem.evaluate(accepted,coefficient,details=True)
        alternatives.append(dict(theta=accepted.tolist(),objective=detail['objective'],
                                 success=bool(refined.success and np.array_equal(accepted,refined.x)),
                                 message=str(refined.message)))
        history.append(dict(seconds=time.monotonic()-start,objective=detail['objective'],
                            data_score=detail['data_score'],penalty=detail['penalty']))
    alternatives.sort(key=lambda item:item['objective'])
    if not alternatives:
        raise RuntimeError("No finite global disk solution was found")
    theta = np.asarray(alternatives[0]['theta'])
    detail = problem.evaluate(theta,coefficient,details=True)
    flags = []
    if not alternatives[0]['success']:
        flags.append('LOCAL_OPTIMIZER_NOT_CONVERGED')
    if not explored.success:
        flags.append('GLOBAL_SEARCH_BUDGET_REACHED')
    near = np.minimum(theta-problem.bounds[:,0],problem.bounds[:,1]-theta) <= 1e-3*span
    if np.any(near[5:]):
        flags.append('ROTATION_SPEED_BOUNDARY')
    if np.any(problem.bounds[5:,1]-theta[5:] <= 1e-3*span[5:]):
        flags.append('ROTATION_SPEED_CEILING')
    if np.any(near[:5]):
        flags.append('GEOMETRY_BOUNDARY')
    if np.any(detail['radial_padding_fraction'] > 1e-8):
        flags.append('OUTER_CONSTANT_SPEED_PADDING_USED')
    if problem.nnode < 3 and coefficient > 0:
        flags.append('TOO_FEW_NODES_FOR_CURVATURE_PENALTY')
    return dict(theta=theta.tolist(),nodes=problem.nodes.tolist(),coefficient=float(coefficient),
                **{key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in detail.items()},
                bin_ids=problem.bin_ids[problem.selected].tolist(),flags=flags,
                alternatives=alternatives,history=history,rejected_evaluations=rejected,
                elapsed_seconds=time.monotonic()-start,
                global_message=str(explored.message),parameter_order=['east0','north0','PA','inclination','Vsys','U_A...','U_B...'])
