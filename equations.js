// equations.js
// JavaScript port of equations.cpp: a leapfrog integrator for a spaceship
// moving near one fixed celestial object.
// Units are dimensionless. The object's strength is the single parameter GM.

const SOFTENING = 0.0001; // same softening length as the C++ version
const SWEEP_DT = 0.005;   // time step used for the impact-parameter sweep

// Acceleration of the spaceship at (x, y) due to the object at (xStar, yStar)
function gravAccel(x, y, xStar, yStar, GM) {
  const dx = x - xStar;
  const dy = y - yStar;
  const r2 = dx * dx + dy * dy + SOFTENING * SOFTENING;
  const r3 = r2 * Math.sqrt(r2);
  const fac = -GM / r3;
  return { ax: fac * dx, ay: fac * dy };
}

// One leapfrog (kick-drift-kick) step. s = { x, y, vx, vy } is updated in place.
function leapfrogStep(s, xStar, yStar, GM, dt) {
  let a = gravAccel(s.x, s.y, xStar, yStar, GM);
  s.vx += 0.5 * dt * a.ax;
  s.vy += 0.5 * dt * a.ay;

  s.x += dt * s.vx;
  s.y += dt * s.vy;

  a = gravAccel(s.x, s.y, xStar, yStar, GM);
  s.vx += 0.5 * dt * a.ax;
  s.vy += 0.5 * dt * a.ay;
}

// Energy per unit mass. Distance is measured from the celestial object
// (the C++ version measured it from the origin, which is only right when
// the object sits at (0, 0)).
function totalEnergy(s, xStar, yStar, GM) {
  const v2 = s.vx * s.vx + s.vy * s.vy;
  const dist = Math.hypot(s.x - xStar, s.y - yStar);
  return 0.5 * v2 - GM / dist;
}

// Angle in degrees between the velocity before and after the encounter
function deflectionAngle(vx0, vy0, vx1, vy1) {
  const m0 = Math.hypot(vx0, vy0);
  const m1 = Math.hypot(vx1, vy1);
  if (m0 === 0 || m1 === 0) return NaN;
  let c = (vx0 * vx1 + vy0 * vy1) / (m0 * m1);
  c = Math.max(-1, Math.min(1, c)); // guard against rounding just outside [-1, 1]
  return Math.acos(c) * 180 / Math.PI;
}

// Run one full flyby.
// p = { x0, y0, vx0, vy0, xStar, yStar, GM, tEnd, nSteps }
function runSimulation(p) {
  const n = p.nSteps;
  const dt = p.tEnd / n;

  const t = new Float64Array(n + 1);
  const x = new Float64Array(n + 1);
  const y = new Float64Array(n + 1);
  const speed = new Float64Array(n + 1);
  const energyErr = new Float64Array(n + 1);

  const s = { x: p.x0, y: p.y0, vx: p.vx0, vy: p.vy0 };
  const E0 = totalEnergy(s, p.xStar, p.yStar, p.GM);
  const denom = Math.abs(E0) || 1;

  let rMin = Math.hypot(s.x - p.xStar, s.y - p.yStar);
  let maxAbsErr = 0;

  t[0] = 0;
  x[0] = s.x;
  y[0] = s.y;
  speed[0] = Math.hypot(s.vx, s.vy);
  energyErr[0] = 0;

  for (let i = 1; i <= n; i++) {
    leapfrogStep(s, p.xStar, p.yStar, p.GM, dt);
    t[i] = i * dt;
    x[i] = s.x;
    y[i] = s.y;
    speed[i] = Math.hypot(s.vx, s.vy);

    const E = totalEnergy(s, p.xStar, p.yStar, p.GM);
    energyErr[i] = (E - E0) / denom;
    maxAbsErr = Math.max(maxAbsErr, Math.abs(energyErr[i]));

    const r = Math.hypot(s.x - p.xStar, s.y - p.yStar);
    if (r < rMin) rMin = r;
  }

  return {
    t, x, y, speed, energyErr,
    rMin,
    maxAbsErr,
    final: { vx: s.vx, vy: s.vy },
    deflectionDeg: deflectionAngle(p.vx0, p.vy0, s.vx, s.vy)
  };
}

// Sweep the impact parameter b.
// Each run starts at (xStar - startDist, yStar + b) moving in the +x direction.
// q = { xStar, yStar, GM, startDist, speed, bMin, bMax, count }
function sweepImpactParameter(q) {
  const tEnd = 3 * q.startDist / q.speed; // long enough to fly past and well beyond
  const nSteps = Math.max(10, Math.round(tEnd / SWEEP_DT));
  const rows = [];

  for (let k = 0; k < q.count; k++) {
    const b = q.count === 1 ? q.bMin : q.bMin + (q.bMax - q.bMin) * k / (q.count - 1);
    const r = runSimulation({
      x0: q.xStar - q.startDist,
      y0: q.yStar + b,
      vx0: q.speed,
      vy0: 0,
      xStar: q.xStar,
      yStar: q.yStar,
      GM: q.GM,
      tEnd,
      nSteps
    });
    rows.push({
      b,
      rMin: r.rMin,
      thetaDeg: r.deflectionDeg,
      vInitial: r.speed[0],
      vFinal: r.speed[nSteps]
    });
  }
  return rows;
}

// Lets the same file be loaded by Node for testing; ignored in the browser.
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    gravAccel, leapfrogStep, totalEnergy, deflectionAngle,
    runSimulation, sweepImpactParameter
  };
}