'use strict';
/* Satellite Measurement Utility — free area and distance measuring on satellite imagery.
   No account required. Everything lives in this file. */

/* ------------------------------------------------------------------ *
 * 1. Geodesy
 *
 * Every measurement is done in a local East/North metre frame anchored
 * at the shape's own centre, built from the true WGS84 radii of
 * curvature at that latitude. Plain planar geometry in that frame is
 * accurate to far better than 0.01% at building scale — tighter than
 * the spherical-earth formula most web mapping tools use, which runs
 * ~0.1% high at mid latitudes.
 * ------------------------------------------------------------------ */
const D2R = Math.PI / 180;
const WGS84_A = 6378137.0;
const WGS84_E2 = (1 / 298.257223563) * (2 - 1 / 298.257223563);

function wrapLon(d) { while (d > 180) d -= 360; while (d < -180) d += 360; return d; }

function frameAt(lat0, lon0) {
  const s = Math.sin(lat0 * D2R), c = Math.cos(lat0 * D2R);
  const w = 1 - WGS84_E2 * s * s, sw = Math.sqrt(w);
  const N = WGS84_A / sw;                        // prime vertical radius
  const M = WGS84_A * (1 - WGS84_E2) / (w * sw); // meridional radius
  const kx = N * c * D2R, ky = M * D2R;          // metres per degree
  return {
    toXY: (lat, lng) => [wrapLon(lng - lon0) * kx, (lat - lat0) * ky],
    toLL: (x, y) => [lat0 + y / ky, lon0 + x / kx]
  };
}

function centreOf(pts) {
  let a = 0, b = 0;
  for (const p of pts) { a += p[0]; b += p[1]; }
  return [a / pts.length, b / pts.length];
}

function pointSegDist(p, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1];
  const L2 = dx * dx + dy * dy;
  let t = L2 > 0 ? ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2 : 0;
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  return Math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
}

function distToPath(p, segs) {
  let m = Infinity;
  for (const s of segs) { const d = pointSegDist(p, s[0], s[1]); if (d < m) m = d; }
  return m;
}

/** Andrew monotone chain convex hull. */
function hull(pts) {
  if (pts.length < 3) return pts.slice();
  const p = pts.slice().sort((u, v) => u[0] - v[0] || u[1] - v[1]);
  const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const build = (src) => {
    const st = [];
    for (const q of src) {
      while (st.length >= 2 && cross(st[st.length - 2], st[st.length - 1], q) <= 0) st.pop();
      st.push(q);
    }
    st.pop();
    return st;
  };
  return build(p).concat(build(p.reverse()));
}

/** Smallest-area enclosing rectangle (rotating calipers over hull edges). */
function minRect(xy) {
  const h = hull(xy);
  if (h.length < 3) return null;
  let best = null;
  for (let i = 0; i < h.length; i++) {
    const a = h[i], b = h[(i + 1) % h.length];
    let ux = b[0] - a[0], uy = b[1] - a[1];
    const L = Math.hypot(ux, uy);
    if (L < 1e-9) continue;
    ux /= L; uy /= L;
    let t0 = Infinity, t1 = -Infinity, s0 = Infinity, s1 = -Infinity;
    for (const q of h) {
      const t = q[0] * ux + q[1] * uy;
      const s = -q[0] * uy + q[1] * ux;
      if (t < t0) t0 = t; if (t > t1) t1 = t;
      if (s < s0) s0 = s; if (s > s1) s1 = s;
    }
    const w = t1 - t0, dd = s1 - s0, area = w * dd;
    if (!best || area < best.area) best = { area, long: Math.max(w, dd), short: Math.min(w, dd) };
  }
  return best;
}

/** Shoelace sums for one closed ring, in a frame already chosen. */
function ringStats(xy) {
  const n = xy.length;
  let twice = 0, perim = 0, cx = 0, cy = 0;
  for (let i = 0; i < n; i++) {
    const a = xy[i], b = xy[(i + 1) % n];
    const cr = a[0] * b[1] - b[0] * a[1];
    twice += cr;
    cx += (a[0] + b[0]) * cr;
    cy += (a[1] + b[1]) * cr;
    perim += Math.hypot(b[0] - a[0], b[1] - a[1]);
  }
  return { twice, perim, cx, cy };
}

/* A ring's centroid falls in its own hole, which would put the label on top of
 * whatever the ring was measured around. This puts it in the band instead: the
 * middle of the longest outer edge, stepped inward half way to the hole. */
function bandAnchor(outer, holes) {
  const segs = [];
  for (const h of holes) for (let i = 0; i < h.length; i++) segs.push([h[i], h[(i + 1) % h.length]]);
  if (!segs.length) return null;
  let best = 0, bestLen = -1;
  for (let i = 0; i < outer.length; i++) {
    const a = outer[i], b = outer[(i + 1) % outer.length];
    const L = Math.hypot(b[0] - a[0], b[1] - a[1]);
    if (L > bestLen) { bestLen = L; best = i; }
  }
  if (bestLen < 1e-9) return null;
  const a = outer[best], b = outer[(best + 1) % outer.length];
  const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
  const gap = distToPath(m, segs);
  if (!(gap > 0)) return m;
  // whichever normal walks towards the hole is the one pointing into the band
  let nx = -(b[1] - a[1]) / bestLen, ny = (b[0] - a[0]) / bestLen;
  if (distToPath([m[0] - nx * gap, m[1] - ny * gap], segs) <
      distToPath([m[0] + nx * gap, m[1] + ny * gap], segs)) { nx = -nx; ny = -ny; }
  return [m[0] + nx * gap / 2, m[1] + ny * gap / 2];
}

/** Area shape: area (m^2), perimeter (m), bounding dimensions, label anchor.
    `holes` are rings cut out of it -- the ground an offset strip does not
    cover. They come off the area and add their own edge to the perimeter. */
function measureArea(pts, holes) {
  const n = pts.length;
  if (n < 2) return { area: 0, perim: 0, dims: null, anchor: pts[0] || [0, 0] };
  const c = centreOf(pts);
  const f = frameAt(c[0], c[1]);
  const xy = pts.map(p => f.toXY(p[0], p[1]));
  const s = ringStats(xy);
  let area = Math.abs(s.twice) / 2, perim = s.perim;

  const holeXY = (holes || []).filter(h => h && h.length >= 3)
                              .map(h => h.map(p => f.toXY(p[0], p[1])));
  for (const h of holeXY) {
    const hs = ringStats(h);
    area -= Math.abs(hs.twice) / 2;
    perim += hs.perim;
  }
  if (area < 0) area = 0;

  let anchor = Math.abs(s.twice) > 1e-9 ? f.toLL(s.cx / (3 * s.twice), s.cy / (3 * s.twice)) : c;
  if (holeXY.length) {
    const a = bandAnchor(xy, holeXY);
    if (a) anchor = f.toLL(a[0], a[1]);
  }
  return { area, perim, dims: n >= 3 ? minRect(xy) : null, anchor };
}

/** Line shape: total length (m), per-segment lengths, label anchor at midpoint. */
function measureLine(pts) {
  const n = pts.length;
  if (n < 2) return { len: 0, segs: [], anchor: pts[0] || [0, 0] };
  const c = centreOf(pts);
  const f = frameAt(c[0], c[1]);
  const xy = pts.map(p => f.toXY(p[0], p[1]));
  const segs = [];
  let len = 0;
  for (let i = 0; i < n - 1; i++) {
    const d = Math.hypot(xy[i + 1][0] - xy[i][0], xy[i + 1][1] - xy[i][1]);
    segs.push(d);
    len += d;
  }
  // anchor halfway along the run
  let acc = 0, anchor = pts[0];
  for (let i = 0; i < segs.length; i++) {
    if (acc + segs[i] >= len / 2) {
      const t = segs[i] > 0 ? (len / 2 - acc) / segs[i] : 0;
      anchor = f.toLL(xy[i][0] + (xy[i + 1][0] - xy[i][0]) * t,
                      xy[i][1] + (xy[i + 1][1] - xy[i][1]) * t);
      break;
    }
    acc += segs[i];
  }
  return { len, segs, anchor };
}

const measureOf = (pts, kind, holes) =>
  kind === 'line' ? measureLine(pts) : measureArea(pts, holes);

/* ------------------------------------------------------------------ *
 * 1b. Outward offset (buffer)
 *
 * Offsets every edge outward, rounds the convex corners with arcs and
 * miters the reflex ones, then discards any candidate point that ended
 * up closer to the original shape than the offset distance. That last
 * filter is what removes the self-intersecting loops a naive offset
 * produces at concave corners.
 * ------------------------------------------------------------------ */
function offsetGeometry(pts, kind, d) {
  if (!(d > 0) || pts.length < 2) return null;
  const c = centreOf(pts);
  const f = frameAt(c[0], c[1]);
  let xy = pts.map(p => f.toXY(p[0], p[1]));

  let ring;
  if (kind === 'line') {
    // walk out and back, so the right of travel is the outside on both passes
    ring = xy.concat(xy.slice(1, -1).reverse());
  } else {
    let s = 0;
    for (let i = 0; i < xy.length; i++) {
      const a = xy[i], b = xy[(i + 1) % xy.length];
      s += a[0] * b[1] - b[0] * a[1];
    }
    ring = s < 0 ? xy.slice().reverse() : xy;   // make it counter-clockwise
  }

  const n = ring.length;
  const segs = [];
  for (let i = 0; i < n; i++) segs.push([ring[i], ring[(i + 1) % n]]);

  const ARC = Math.PI / 32;   // arc step; finer than this is invisible on screen
  const cand = [];
  for (let i = 0; i < n; i++) {
    const a = ring[i], b = ring[(i + 1) % n];
    const dx = b[0] - a[0], dy = b[1] - a[1];
    const L = Math.hypot(dx, dy);
    if (L < 1e-9) continue;
    const nx = dy / L, ny = -dx / L;             // outward normal (right of travel)
    cand.push([a[0] + nx * d, a[1] + ny * d]);
    cand.push([b[0] + nx * d, b[1] + ny * d]);

    // direction of the next non-degenerate edge, to know how to turn the corner
    let nextPt = null;
    for (let k = 1; k <= n; k++) {
      const q = ring[(i + 1 + k) % n];
      if (Math.hypot(q[0] - b[0], q[1] - b[1]) > 1e-9) { nextPt = q; break; }
    }
    if (!nextPt) continue;
    const ex = nextPt[0] - b[0], ey = nextPt[1] - b[1];
    const L2 = Math.hypot(ex, ey);
    const mx = ey / L2, my = -ex / L2;

    const a0 = Math.atan2(ny, nx);
    let sweep = Math.atan2(my, mx) - a0;
    while (sweep <= -Math.PI) sweep += 2 * Math.PI;
    while (sweep > Math.PI) sweep -= 2 * Math.PI;

    if (sweep > 1e-9) {                          // convex — round the corner
      const steps = Math.max(1, Math.ceil(sweep / ARC));
      for (let s = 1; s < steps; s++) {
        const t = a0 + sweep * s / steps;
        cand.push([b[0] + Math.cos(t) * d, b[1] + Math.sin(t) * d]);
      }
    } else if (sweep < -1e-9) {                  // reflex — where the offsets cross
      const half = sweep / 2;
      if (Math.abs(half) < Math.PI * 0.49) {
        const r = d / Math.cos(half);
        if (r < d * 6) {
          const t = a0 + half;
          cand.push([b[0] + Math.cos(t) * r, b[1] + Math.sin(t) * r]);
        }
      }
    }
  }

  const tol = d * 1e-3;
  return closeRing(cand.filter(p => distToPath(p, segs) >= d - tol), tol, f);
}

/* Drops the repeats an offset leaves behind and hands the ring back as
   lat/lngs, unclosed, the way every shape in this app is held. */
function closeRing(xy, tol, f) {
  const clean = [];
  for (const q of xy) {
    const last = clean[clean.length - 1];
    if (!last || Math.hypot(q[0] - last[0], q[1] - last[1]) > tol) clean.push(q);
  }
  if (clean.length > 2) {
    const a = clean[0], b = clean[clean.length - 1];
    if (Math.hypot(a[0] - b[0], a[1] - b[1]) <= tol) clean.pop();
  }
  if (clean.length < 3) return null;
  return clean.map(q => f.toLL(q[0], q[1]));
}

/* Offsets only the sides named in `sides`, which are edge indices into the
 * original point list: edge i runs from pts[i] to pts[i+1].
 *
 * The full buffer above rounds its corners because every edge moves by the
 * same distance and the corner has to be filled in. Here neighbouring edges
 * can move by different amounts, so the corner is simply where the two moved
 * edge lines now cross -- a mitre. That is also what a setback on two sides of
 * a building looks like on the ground: the walls stay straight and meet at a
 * point. A corner between edges that both stayed put does not move at all.
 *
 * A very acute corner throws the mitre a long way out, so past four times the
 * offset it is cut off square instead.
 */
function offsetSides(pts, d, sides) {
  const r = sideRing(pts, d, sides);
  if (!r) return null;
  const out = [];
  for (const corner of r.corners) for (const q of corner) out.push(q);
  return closeRing(out, Math.max(d * 1e-6, 1e-9), r.f);
}

/* The shape re-expressed counter-clockwise, so the right of travel is always
   the outside, together with the mitred corners for this set of sides. On a
   reversed ring, ring edge i is the original edge n-2-i -- side i has to keep
   pointing at the edge from pts[i] to pts[i+1] whichever way the ring was
   drawn, because that is what the numbered chips in the sidebar name. */
function sideRing(pts, d, sides) {
  if (!(d > 0) || pts.length < 3 || !sides || !sides.size) return null;
  const c = centreOf(pts);
  const f = frameAt(c[0], c[1]);
  const xy = pts.map(p => f.toXY(p[0], p[1]));
  const n = xy.length;

  let s = 0;
  for (let i = 0; i < n; i++) {
    const a = xy[i], b = xy[(i + 1) % n];
    s += a[0] * b[1] - b[0] * a[1];
  }
  const ccw = s >= 0;
  const ring = ccw ? xy : xy.slice().reverse();
  const srcOf = (i) => ccw ? i : (n - 2 - i + n) % n;
  const moved = (i) => sides.has(srcOf(i));

  const corners = mitredCorners(ring, d, moved);
  return corners && { f, ring, n, moved, corners };
}

/* Where each ring vertex ends up once the moved edges have moved: the crossing
   of the two edge lines either side of it, one or both of them displaced. A
   corner between two edges that both stayed put does not move at all. A very
   acute corner throws the mitre a long way out, so past four times the offset
   it is cut off square, which is the one case where a vertex gives two points
   instead of one.
   Kept per vertex so a caller can walk the offset alongside the original,
   which is what offsetStrip does. */
function mitredCorners(ring, d, moved) {
  const n = ring.length;
  const lines = [];
  for (let i = 0; i < n; i++) {
    const a = ring[i], b = ring[(i + 1) % n];
    const dx = b[0] - a[0], dy = b[1] - a[1];
    const L = Math.hypot(dx, dy);
    if (L < 1e-9) return null;
    const nx = dy / L, ny = -dx / L;
    const push = moved(i) ? d : 0;
    lines.push({
      v: [dx / L, dy / L],
      a: [a[0] + nx * push, a[1] + ny * push],
      b: [b[0] + nx * push, b[1] + ny * push]
    });
  }

  const corners = [];
  for (let i = 0; i < n; i++) {
    const prev = lines[(i - 1 + n) % n], cur = lines[i];
    const cross = prev.v[0] * cur.v[1] - prev.v[1] * cur.v[0];
    if (Math.abs(cross) < 1e-9) { corners.push([cur.a]); continue; }   // straight through
    const wx = cur.a[0] - prev.a[0], wy = cur.a[1] - prev.a[1];
    const t = (wx * cur.v[1] - wy * cur.v[0]) / cross;
    const hit = [prev.a[0] + prev.v[0] * t, prev.a[1] + prev.v[1] * t];
    corners.push(Math.hypot(hit[0] - ring[i][0], hit[1] - ring[i][1]) > d * 4
      ? [prev.b, cur.a]                                                // bevel
      : [hit]);
  }
  return corners;
}

/* The same outward offset as offsetGeometry, but with every corner mitred
 * instead of rounded: the sides stay straight and meet at a point.
 *
 * The rounded buffer is the more faithful answer -- it is the set of points
 * within d of the shape -- but it spends about ninety vertices on the corners
 * of a warehouse, and a shape with ninety corners cannot be nudged by hand.
 * A mitre gives back a ring with as many corners as the building has, which is
 * also what a setback actually looks like on the ground: the walls stay
 * straight and meet at a point.
 *
 * A notch narrower than twice the offset folds its corners back inside the
 * shape. Those are dropped by the same test offsetGeometry uses -- a corner
 * that ended up closer to the original than the offset distance was never on
 * the offset in the first place.
 */
function offsetMitred(pts, d) {
  const r = sideRing(pts, d, new Set(pts.map((_, i) => i)));
  if (!r) return null;
  const { f, ring, corners } = r;
  const segs = ring.map((p, i) => [p, ring[(i + 1) % ring.length]]);
  const tol = d * 1e-3;

  const all = [];
  for (const corner of corners) for (const q of corner) all.push(q);
  const keep = all.filter(q => distToPath(q, segs) >= d - tol);
  return closeRing(keep.length >= 3 ? keep : all, tol, f);
}

/* Just the ground the offset adds along the chosen sides, as a shape of its
 * own: the moved edges on the outside, the original edges on the inside. The
 * building it was measured from is left alone, so nothing is counted twice.
 *
 * Sides that touch make one band -- pushing two adjacent walls out gives a
 * single L. Sides on opposite walls make two, because that is what they are on
 * the ground, so this hands back a list of rings rather than one.
 *
 * A run of moved edges i0..i1 spans vertices i0..i1+1. Walking those vertices
 * forward along the original and then back across their mitred positions
 * closes the band exactly, with no polygon clipping involved: each end of the
 * run is the corner where a moved edge meets one that stayed where it was.
 */
function offsetStrip(pts, d, sides) {
  const r = sideRing(pts, d, sides);
  if (!r) return null;
  const { f, ring, n, moved, corners } = r;

  // Start at the first moved edge whose predecessor stayed put, so that no run
  // wraps past the join. With every edge moved there is no such edge and the
  // band is a closed ring -- that is offsetGeometry's job, not this one.
  let start = -1;
  for (let i = 0; i < n; i++) if (moved(i) && !moved((i - 1 + n) % n)) { start = i; break; }
  if (start < 0) return null;

  const runs = [];
  let run = null;
  for (let k = 0; k < n; k++) {
    const i = (start + k) % n;
    if (!moved(i)) { run = null; continue; }
    if (!run) { run = []; runs.push(run); }
    run.push(i);
  }

  const tol = Math.max(d * 1e-6, 1e-9);
  const out = [];
  for (const edges of runs) {
    // Counted forward from the run's first vertex, so that a run which wraps
    // past vertex 0 still walks in one direction.
    const i0 = edges[0], end = i0 + edges.length;
    const band = [];
    for (let v = i0; v <= end; v++) band.push(ring[v % n]);              // inside, forward
    for (let v = end; v >= i0; v--) {                                    // outside, back
      const corner = corners[v % n];
      for (let k = corner.length - 1; k >= 0; k--) band.push(corner[k]);
    }
    const llr = closeRing(band, tol, f);
    if (llr) out.push(llr);
  }
  return out.length ? out : null;
}

/* The length of each side, in metres, in the order the sides are numbered. */
function sideLengths(pts) {
  if (pts.length < 2) return [];
  const c = centreOf(pts);
  const f = frameAt(c[0], c[1]);
  const xy = pts.map(p => f.toXY(p[0], p[1]));
  const out = [];
  for (let i = 0; i < xy.length; i++) {
    const a = xy[i], b = xy[(i + 1) % xy.length];
    out.push(Math.hypot(b[0] - a[0], b[1] - a[1]));
  }
  return out;
}

/* ------------------------------------------------------------------ *
 * 1c. Union
 *
 * Merging measured shapes into one is a polygon union. This does it the
 * direct way rather than with a clipping library, because the app carries no
 * dependencies and the shapes here are small:
 *
 *   1. cut every edge at every crossing with every other edge, so that no two
 *      edges meet anywhere but at a shared end;
 *   2. throw away the pieces that are not on the outside of the union -- a
 *      piece with the union on both sides of it, or on neither, is not part
 *      of its boundary. Which side is which is decided by sampling a point
 *      just off each side, so a shared edge drops out whichever direction the
 *      two shapes were drawn in;
 *   3. point what is left so the union is always on its left, and walk the
 *      loops out of it.
 *
 * The loops that come back anticlockwise are outlines and the clockwise ones
 * are holes, which is exactly the shape model everything else here uses.
 *
 * All of it happens in one local metre frame shared by every input, so it is
 * plain planar geometry and the result measures with the same shoelace as
 * everything else.
 * ------------------------------------------------------------------ */
const UNION_TOL = 1e-6;      // metres; two corners this close are one corner

/** Ray casting. Orientation does not matter. */
function pointInRing(p, ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const a = ring[i], b = ring[j];
    if ((a[1] > p[1]) !== (b[1] > p[1]) &&
        p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) inside = !inside;
  }
  return inside;
}

/** Inside the outline and outside every hole. */
function pointInPoly(p, poly) {
  if (!pointInRing(p, poly.outer)) return false;
  for (const h of poly.holes) if (pointInRing(p, h)) return false;
  return true;
}

/* Where segment s has to be cut so that it meets t only at an end. Proper
   crossings give one parameter on each; two collinear segments that overlap
   give each other's ends. */
function addCuts(s, t, cs, ct) {
  const rx = s[1][0] - s[0][0], ry = s[1][1] - s[0][1];
  const ux = t[1][0] - t[0][0], uy = t[1][1] - t[0][1];
  const ls = Math.hypot(rx, ry), lt = Math.hypot(ux, uy);
  if (ls < UNION_TOL || lt < UNION_TOL) return;
  const wx = t[0][0] - s[0][0], wy = t[0][1] - s[0][1];
  const den = rx * uy - ry * ux;

  if (Math.abs(den) > 1e-12 * ls * lt) {
    const a = (wx * uy - wy * ux) / den;
    const b = (wx * ry - wy * rx) / den;
    if (a > -UNION_TOL / ls && a < 1 + UNION_TOL / ls &&
        b > -UNION_TOL / lt && b < 1 + UNION_TOL / lt) {
      cs.push(Math.min(1, Math.max(0, a)));
      ct.push(Math.min(1, Math.max(0, b)));
    }
    return;
  }
  if (Math.abs(wx * ry - wy * rx) / ls > UNION_TOL) return;   // parallel, not collinear
  for (const q of t) {
    const a = ((q[0] - s[0][0]) * rx + (q[1] - s[0][1]) * ry) / (ls * ls);
    if (a > 0 && a < 1) cs.push(a);
  }
  for (const q of s) {
    const b = ((q[0] - t[0][0]) * ux + (q[1] - t[0][1]) * uy) / (lt * lt);
    if (b > 0 && b < 1) ct.push(b);
  }
}

/* Corners within UNION_TOL of one another are the same corner. Snapping to the
   nearest existing one rather than to a grid keeps two points that straddle a
   grid line from being pulled apart. */
function nodeIndex(tol) {
  const cell = tol * 4;
  const buckets = new Map();
  const pts = [];
  return {
    pts,
    id(p) {
      const cx = Math.floor(p[0] / cell), cy = Math.floor(p[1] / cell);
      for (let dx = -1; dx <= 1; dx++) {
        for (let dy = -1; dy <= 1; dy++) {
          const b = buckets.get((cx + dx) + ':' + (cy + dy));
          if (!b) continue;
          for (const i of b) {
            if (Math.hypot(pts[i][0] - p[0], pts[i][1] - p[1]) <= tol) return i;
          }
        }
      }
      const i = pts.length;
      pts.push(p);
      const k = cx + ':' + cy;
      if (!buckets.has(k)) buckets.set(k, []);
      buckets.get(k).push(i);
      return i;
    }
  };
}

/* Merges shapes into as few as they make. Takes and returns {pts, holes} in
   lat/lng, holes optional. Returns null if there is nothing usable in them. */
function unionShapes(shapes) {
  const rings = [];
  for (const sh of shapes) {
    if (!sh.pts || sh.pts.length < 3) continue;
    rings.push(sh.pts, ...(sh.holes || []).filter(h => h && h.length >= 3));
  }
  if (!rings.length) return null;
  const c = centreOf(rings.flat());
  const f = frameAt(c[0], c[1]);
  const toXY = (r) => r.map(q => f.toXY(q[0], q[1]));

  const polys = shapes
    .filter(sh => sh.pts && sh.pts.length >= 3)
    .map(sh => ({ outer: toXY(sh.pts),
                  holes: (sh.holes || []).filter(h => h && h.length >= 3).map(toXY) }));
  const inUnion = (q) => polys.some(poly => pointInPoly(q, poly));

  const segs = [];
  for (const poly of polys) {
    for (const ring of [poly.outer, ...poly.holes]) {
      for (let i = 0; i < ring.length; i++) segs.push([ring[i], ring[(i + 1) % ring.length]]);
    }
  }

  const cuts = segs.map(() => [0, 1]);
  for (let i = 0; i < segs.length; i++) {
    for (let j = i + 1; j < segs.length; j++) addCuts(segs[i], segs[j], cuts[i], cuts[j]);
  }

  const idx = nodeIndex(UNION_TOL);
  const seen = new Set();
  const pieces = [];
  for (let i = 0; i < segs.length; i++) {
    const [p, q] = segs[i];
    const ts = cuts[i].filter(v => v >= 0 && v <= 1).sort((a, b) => a - b);
    const at = (t) => [p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t];
    for (let k = 0; k + 1 < ts.length; k++) {
      const a = idx.id(at(ts[k])), b = idx.id(at(ts[k + 1]));
      if (a === b) continue;
      // two shapes sharing a wall contribute the same piece twice; one is enough
      const key = Math.min(a, b) + ':' + Math.max(a, b);
      if (seen.has(key)) continue;
      seen.add(key);
      pieces.push([a, b]);
    }
  }

  const edges = [];
  for (const [a, b] of pieces) {
    const p = idx.pts[a], q = idx.pts[b];
    const dx = q[0] - p[0], dy = q[1] - p[1];
    const len = Math.hypot(dx, dy);
    if (len < UNION_TOL) continue;
    const e = Math.min(len / 4, 0.01);
    const mx = (p[0] + q[0]) / 2, my = (p[1] + q[1]) / 2;
    const nx = -dy / len * e, ny = dx / len * e;          // left of a -> b
    const left = inUnion([mx + nx, my + ny]);
    const right = inUnion([mx - nx, my - ny]);
    if (left === right) continue;                          // inside or outside on both sides
    edges.push(left ? { a, b, dx, dy } : { a: b, b: a, dx: -dx, dy: -dy });
  }
  if (!edges.length) return null;

  const outAt = new Map();
  edges.forEach((e, i) => {
    if (!outAt.has(e.a)) outAt.set(e.a, []);
    outAt.get(e.a).push(i);
  });

  /* Where two boundaries meet at a single corner there is a choice of exits.
     Taking the first one clockwise of the way we came keeps the traversal
     hugging the same side of the boundary all the way round. */
  const nextFrom = (i, used) => {
    const back = Math.atan2(-edges[i].dy, -edges[i].dx);
    let best = -1, bestAng = Infinity;
    for (const k of outAt.get(edges[i].b) || []) {
      if (used[k]) continue;
      let ang = back - Math.atan2(edges[k].dy, edges[k].dx);
      while (ang <= 1e-12) ang += 2 * Math.PI;
      while (ang > 2 * Math.PI) ang -= 2 * Math.PI;
      if (ang < bestAng) { bestAng = ang; best = k; }
    }
    return best;
  };

  const used = new Array(edges.length).fill(false);
  const loops = [];
  for (let start = 0; start < edges.length; start++) {
    if (used[start]) continue;
    const loop = [];
    for (let e = start; e >= 0 && !used[e]; e = nextFrom(e, used)) {
      used[e] = true;
      loop.push(idx.pts[edges[e].a]);
    }
    if (loop.length >= 3) loops.push({ ring: loop, twice: ringStats(loop).twice });
  }

  const outers = loops.filter(l => l.twice > 1e-9).map(l => ({ ring: l.ring, area: l.twice / 2, holes: [] }));
  const holes = loops.filter(l => l.twice < -1e-9);
  if (!outers.length) return null;
  for (const h of holes) {
    // an edge midpoint, not a corner: a corner may be shared with the outline
    const m = [(h.ring[0][0] + h.ring[1][0]) / 2, (h.ring[0][1] + h.ring[1][1]) / 2];
    let host = null;
    for (const o of outers) {
      if (pointInRing(m, o.ring) && (!host || o.area < host.area)) host = o;
    }
    if (host) host.holes.push(h.ring);
  }

  const back = (r) => dropCollinear(r).map(q => f.toLL(q[0], q[1]));
  return outers
    .sort((a, b) => b.area - a.area)
    .map(o => ({ pts: back(o.ring), holes: o.holes.length ? o.holes.map(back) : null }));
}

/* Two shapes that shared a wall leave the ends of that wall behind as corners
   sitting in the middle of a straight side. They measure the same either way,
   but the whole reason for merging is to end up with a shape you can still
   take hold of, so they go. */
function dropCollinear(ring) {
  const out = [];
  for (let i = 0; i < ring.length; i++) {
    const a = ring[(i - 1 + ring.length) % ring.length], b = ring[i], c = ring[(i + 1) % ring.length];
    const ux = b[0] - a[0], uy = b[1] - a[1], vx = c[0] - b[0], vy = c[1] - b[1];
    const len = Math.hypot(vx, vy);
    if (len > UNION_TOL && Math.abs(ux * vy - uy * vx) / len <= UNION_TOL) continue;
    out.push(b);
  }
  return out.length >= 3 ? out : ring;
}

/* ------------------------------------------------------------------ *
 * 2. Units
 * ------------------------------------------------------------------ */
const AREA_UNITS = {
  ft2:  { label: 'sq ft',   per_m2: 10.763910416709722, dp: 0 },
  m2:   { label: 'sq m',    per_m2: 1,                  dp: 0 },
  acre: { label: 'acres',   per_m2: 1 / 4046.8564224,   dp: 3 },
  ha:   { label: 'ha',      per_m2: 1 / 10000,          dp: 3 },
  yd2:  { label: 'sq yd',   per_m2: 1 / 0.83612736,     dp: 0 },
  mi2:  { label: 'sq mi',   per_m2: 1 / 2589988.110336, dp: 4 }
};
const LEN_UNITS = {
  ft: { label: 'ft', per_m: 3.280839895013123, dp: 1 },
  m:  { label: 'm',  per_m: 1,                 dp: 1 },
  yd: { label: 'yd', per_m: 1 / 0.9144,        dp: 1 }
};

const num = (v, dp) => v.toLocaleString(undefined, { minimumFractionDigits: dp, maximumFractionDigits: dp });
function fmtArea(m2, withUnit = true) {
  const u = AREA_UNITS[state.areaUnit];
  return num(m2 * u.per_m2, u.dp) + (withUnit ? ' ' + u.label : '');
}
function fmtLen(m, withUnit = true) {
  const u = LEN_UNITS[state.lenUnit];
  return num(m * u.per_m, u.dp) + (withUnit ? ' ' + u.label : '');
}
const pitchFactor = () => Math.sqrt(1 + Math.pow(state.pitch / 12, 2));

/* ------------------------------------------------------------------ *
 * 3. State
 * ------------------------------------------------------------------ */
const COLOR = { add: '#38bdf8', subtract: '#f87171', line: '#fbbf24' };
const STORE_KEY = 'satellite-measurement-utility/v1';
const LEGACY_KEY = 'warehouse-measure/v1';   // renamed; keep older saves loadable

const state = {
  shapes: [],        // {id,name,kind:'area'|'line',mode,pts,layer,label,m}
  selected: null,
  areaUnit: 'ft2',
  lenUnit: 'ft',
  pitch: 0,
  rate: null,
  labels: true,
  sideNums: false,
  offsetDist: 100,
  captureMargin: 150,
  gkey: ''
};
let nextId = 1;

const $ = (id) => document.getElementById(id);

/* ------------------------------------------------------------------ *
 * 4. Map and imagery
 * ------------------------------------------------------------------ */
const esri = L.tileLayer(
  'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  { maxZoom: 22, maxNativeZoom: 20, attribution: 'Imagery &copy; Esri, Maxar, Earthstar Geographics' });

const labelsOverlay = L.tileLayer(
  'https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile/{z}/{y}/{x}',
  { maxZoom: 22, maxNativeZoom: 20, attribution: '&copy; Esri' });

const usgs = L.tileLayer.wms(
  'https://basemap.nationalmap.gov/arcgis/services/USGSImageryOnly/MapServer/WMSServer',
  { layers: '0', format: 'image/jpeg', transparent: false, maxZoom: 22,
    attribution: 'USGS National Map (public domain)' });

const osm = L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',
  { maxZoom: 19, attribution: '&copy; OpenStreetMap contributors' });

const map = L.map('map', {
  center: [39.5, -98.35], zoom: 5, maxZoom: 22,
  layers: [esri], zoomControl: true, doubleClickZoom: true, attributionControl: true
});
L.control.scale({ imperial: true, metric: true, position: 'bottomright' }).addTo(map);

// edit handles need to sit above the measurement labels, which live in tooltipPane (650)
map.createPane('handles');
map.getPane('handles').style.zIndex = 680;

const baseLayers = { 'Satellite (Esri)': esri, 'Satellite (USGS, US only)': usgs, 'Street map': osm };
const layerControl = L.control.layers(baseLayers, { 'Place labels': labelsOverlay },
                                      { position: 'topright' }).addTo(map);

/* --- Google Map Tiles API. Optional: needs the user's own key. --------
   Google's imagery may only be served through their own endpoints, so this
   uses the official Map Tiles API — create a session, then pull tiles with
   it. Any failure drops back to Esri rather than leaving a blank map. */
let googleLayer = null;

async function enableGoogle(key) {
  const msg = (t) => { $('gkeyMsg').textContent = t; };
  msg('Contacting Google…');
  let session;
  try {
    const r = await fetch('https://tile.googleapis.com/v1/createSession?key=' + encodeURIComponent(key), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ mapType: 'satellite', language: 'en-US', region: 'US' })
    });
    if (!r.ok) {
      let detail = 'HTTP ' + r.status;
      try { const j = await r.json(); if (j.error && j.error.message) detail = j.error.message; } catch (e) {}
      throw new Error(detail);
    }
    session = (await r.json()).session;
    if (!session) throw new Error('no session token returned');
  } catch (e) {
    msg('Google unavailable, staying on Esri — ' + e.message);
    return false;
  }

  disableGoogle();
  googleLayer = L.tileLayer(
    `https://tile.googleapis.com/v1/2dtiles/{z}/{x}/{y}?session=${encodeURIComponent(session)}&key=${encodeURIComponent(key)}`,
    { maxZoom: 22, maxNativeZoom: 21, attribution: 'Imagery &copy; Google' });

  let errors = 0;
  googleLayer.on('tileerror', () => {
    if (++errors === 8) {
      msg('Google tiles kept failing, switched back to Esri.');
      disableGoogle();
      map.addLayer(esri);
    }
  });

  baseLayers['Satellite (Google)'] = googleLayer;
  layerControl.addBaseLayer(googleLayer, 'Satellite (Google)');
  map.removeLayer(esri);
  map.addLayer(googleLayer);
  msg('Using Google imagery.');
  return true;
}

function disableGoogle() {
  if (!googleLayer) return;
  map.removeLayer(googleLayer);
  layerControl.removeLayer(googleLayer);
  delete baseLayers['Satellite (Google)'];
  googleLayer = null;
}

/* ------------------------------------------------------------------ *
 * 5. Shapes
 * ------------------------------------------------------------------ */
function defaultName(kind) {
  const n = state.shapes.filter(s => s.kind === kind).length + 1;
  return kind === 'line' ? `Line ${n}` : `Building ${n}`;
}

function addShape(pts, opts = {}) {
  const kind = opts.kind === 'line' ? 'line' : 'area';
  const holes = kind === 'area' && opts.holes && opts.holes.length
    ? opts.holes.filter(h => h && h.length >= 3) : null;
  const sh = {
    id: nextId++,
    name: opts.name || defaultName(kind),
    kind,
    mode: opts.mode === 'subtract' ? 'subtract' : 'add',
    color: opts.color || null,
    // Hidden shapes stay in the list and in the saved file. They are off the
    // map and out of the totals, which is the whole point of the toggle.
    hidden: !!opts.hidden,
    pts,
    // Ground inside the outline that the shape does not cover -- what an
    // offset strip leaves for the building it was measured around.
    holes: holes && holes.length ? holes : null,
    layer: null,
    label: null
  };
  sh.layer = kind === 'line' ? L.polyline(pts, styleFor(sh))
                             : L.polygon(latLngsOf(sh), styleFor(sh));
  if (!sh.hidden) sh.layer.addTo(map);
  sh.layer.on('click', (e) => {
    L.DomEvent.stop(e);
    if (tool) return;
    const oe = e.originalEvent;
    if (oe && (oe.ctrlKey || oe.metaKey)) toggleMulti(sh.id); else select(sh.id);
  });
  sh.layer.on('mousedown', (e) => startMove(sh, e));
  state.shapes.push(sh);
  refreshShape(sh);
  return sh;
}

/** A shape's own colour if it has been customised, otherwise the default
    for its kind and mode. */
const colorOf = (sh) => sh.color || (sh.kind === 'line' ? COLOR.line : COLOR[sh.mode]);

/** What Leaflet wants: one ring, or the outline followed by its holes. */
const latLngsOf = (sh) => sh.holes && sh.holes.length ? [sh.pts, ...sh.holes] : sh.pts;

/** Every ring the shape is drawn from, outline first. */
const ringsOf = (sh) => sh.kind === 'line' || !sh.holes ? [sh.pts] : [sh.pts, ...sh.holes];

const measureShape = (sh) => sh.m || measureOf(sh.pts, sh.kind, sh.holes);

/** The shape as it is saved and as undo remembers it -- no Leaflet, no id. */
const plainShape = (sh) => ({ name: sh.name, kind: sh.kind, mode: sh.mode, color: sh.color,
                              hidden: sh.hidden, pts: sh.pts, holes: sh.holes });

/** Slide a whole shape, holes and all. */
function translateShape(sh, dLat, dLng) {
  const shift = (ring) => ring.map(q => [q[0] + dLat, q[1] + dLng]);
  sh.pts = shift(sh.pts);
  if (sh.holes) sh.holes = sh.holes.map(shift);
  sh.layer.setLatLngs(latLngsOf(sh));
}

function styleFor(sh) {
  const c = colorOf(sh);
  const on = state.selected === sh.id || multi.has(sh.id);
  return {
    color: c,
    weight: sh.kind === 'line' ? (on ? 5 : 4) : (on ? 3 : 2),
    opacity: 1,
    fillColor: c,
    fillOpacity: sh.kind === 'line' ? 0 : (on ? 0.3 : 0.18),
    dashArray: sh.mode === 'subtract' && sh.kind !== 'line' ? '6 4' : null
  };
}

function labelFor(sh, m) {
  if (sh.kind === 'line') return `${escapeHtml(sh.name)}<br>${fmtLen(m.len)}`;
  return `${escapeHtml(sh.name)}<br>${sh.mode === 'subtract' ? '&minus;' : ''}${fmtArea(m.area)}`;
}

function refreshShape(sh) {
  sh.layer.setLatLngs(latLngsOf(sh));
  sh.layer.setStyle(styleFor(sh));
  const m = measureOf(sh.pts, sh.kind, sh.holes);
  sh.m = m;
  // Measured either way: a hidden shape still shows its own figure in the
  // list, so it can be judged before being counted again.
  if (sh.hidden && map.hasLayer(sh.layer)) map.removeLayer(sh.layer);
  if (!sh.hidden && !map.hasLayer(sh.layer)) sh.layer.addTo(map);
  if (sh.label) { map.removeLayer(sh.label); sh.label = null; }
  const has = sh.kind === 'line' ? m.len > 0 : m.area > 0;
  if (state.labels && has && !sh.hidden) {
    sh.label = L.tooltip({ permanent: true, direction: 'center', className: 'measure-label', interactive: false })
      .setLatLng(m.anchor).setContent(labelFor(sh, m)).addTo(map);
  }
}

function removeShape(id) {
  const i = state.shapes.findIndex(s => s.id === id);
  if (i < 0) return;
  const sh = state.shapes[i];
  map.removeLayer(sh.layer);
  if (sh.label) map.removeLayer(sh.label);
  state.shapes.splice(i, 1);
  multi.delete(id);
  if (state.selected === id) select(null);
  else { renderAll(); save(); }
}

function select(id) {
  state.selected = id;
  offsetNote = '';
  multi.clear();
  clearHandles();
  state.shapes.forEach(s => s.layer.setStyle(styleFor(s)));
  const sh = state.shapes.find(s => s.id === id);
  // A hidden shape can be selected from the list -- to rename it, or to bring
  // it back -- but it has no outline on the map to hang edit handles off.
  if (sh && !sh.hidden) buildHandles(sh);
  sidePick = { id: sh && !sh.hidden && sh.kind === 'area' ? id : null, sides: null };
  hotSide = null;
  renderAll();
  drawSideHighlight();
  drawSideNums();
}

/** Show or hide one shape. */
function setHidden(sh, hidden) {
  pushUndo();
  sh.hidden = !!hidden;
  if (sh.hidden && state.selected === sh.id) select(null);
  else { refreshShape(sh); renderAll(); save(); }
}

/** Both buttons under the list. */
function setAllHidden(hidden) {
  pushUndo();
  state.shapes.forEach(s => { s.hidden = !!hidden; refreshShape(s); });
  if (hidden) select(null);
  else { renderAll(); save(); }
}

const selectedShape = () => state.shapes.find(s => s.id === state.selected) || null;

/* Merging needs more than one shape, but everything else here -- the handles,
   the offset box, the side chips -- is about exactly one. So there is still a
   selection proper, and Ctrl+click adds others alongside it. */
const multi = new Set();
const pickedShapes = () =>
  state.shapes.filter(s => s.id === state.selected || multi.has(s.id));

function toggleMulti(id) {
  if (state.selected === null || state.selected === undefined) { select(id); return; }
  if (id === state.selected) return;          // the selection proper cannot be dropped
  if (multi.has(id)) multi.delete(id); else multi.add(id);
  state.shapes.forEach(s => s.layer.setStyle(styleFor(s)));
  renderAll();
}

/* ------------------------------------------------------------------ *
 * 6. Vertex handles and dragging a whole shape
 * ------------------------------------------------------------------ */
let handles = [];
function clearHandles() { handles.forEach(h => map.removeLayer(h)); handles = []; }

function handleIcon(mid) {
  const s = mid ? 9 : 12;
  return L.divIcon({ className: 'vtx' + (mid ? ' mid' : ''), iconSize: [s, s], iconAnchor: [s / 2, s / 2] });
}

/** Centre grip. Dragging it slides the whole shape — the only reliable way to
    grab a thin line, which is too narrow a target to hit reliably. */
function addMoveHandle(sh) {
  const mv = L.marker(sh.m.anchor, {
    icon: L.divIcon({ className: 'vtx move', iconSize: [20, 20], iconAnchor: [10, 10], html: '&#10021;' }),
    draggable: true, zIndexOffset: 1100, pane: 'handles', title: 'Drag to move the whole shape'
  }).addTo(map);
  let last = L.latLng(sh.m.anchor);
  mv.on('dragstart', (e) => {
    beginEdit();
    last = e.target.getLatLng();
    handles.forEach(h => { if (h !== mv) map.removeLayer(h); });
    handles = [mv];
  });
  mv.on('drag', (e) => {
    const cur = e.target.getLatLng();
    translateShape(sh, cur.lat - last.lat, cur.lng - last.lng);
    last = cur;
  });
  mv.on('dragend', () => { commitEdit(); refreshShape(sh); renderAll(); save(); buildHandles(sh); });
  handles.push(mv);
}

function buildHandles(sh) {
  clearHandles();
  if (tool) return;
  const closed = sh.kind !== 'line';
  addMoveHandle(sh);
  sh.pts.forEach((p, i) => {
    const mk = L.marker(p, { icon: handleIcon(false), draggable: true, zIndexOffset: 1000,
                             pane: 'handles' }).addTo(map);
    mk.on('dragstart', beginEdit);
    mk.on('drag', (e) => {
      const ll = e.target.getLatLng();
      sh.pts[i] = [ll.lat, ll.lng];
      sh.layer.setLatLngs(latLngsOf(sh));
      liveReadout(measureOf(sh.pts, sh.kind, sh.holes), sh.kind);
    });
    mk.on('dragend', () => {
      commitEdit();
      refreshShape(sh); buildHandles(sh); renderAll(); save(); setReadout('');
    });
    mk.on('contextmenu', (e) => {
      L.DomEvent.stop(e);
      if (sh.pts.length <= (closed ? 3 : 2)) return;
      pushUndo();
      sh.pts.splice(i, 1);
      refreshShape(sh); buildHandles(sh); renderAll(); save();
    });
    handles.push(mk);

    if (!closed && i === sh.pts.length - 1) return;   // open line has no closing midpoint
    const q = sh.pts[(i + 1) % sh.pts.length];
    const mid = [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2];
    const mm = L.marker(mid, { icon: handleIcon(true), zIndexOffset: 900, pane: 'handles',
      title: 'Click to add a corner here' }).addTo(map);
    mm.on('click', (e) => {
      L.DomEvent.stop(e);
      pushUndo();
      sh.pts.splice(i + 1, 0, mid);
      refreshShape(sh); buildHandles(sh); renderAll(); save();
    });
    handles.push(mm);
  });
}

/** Drag anywhere inside a selected shape to slide the whole thing. */
let moving = null;
function startMove(sh, e) {
  if (tool || state.selected !== sh.id) return;
  L.DomEvent.stop(e);
  map.dragging.disable();
  clearHandles();
  beginEdit();
  moving = { sh, last: e.latlng, moved: false };
}
map.on('mousemove', (e) => {
  if (!moving) return;
  translateShape(moving.sh, e.latlng.lat - moving.last.lat, e.latlng.lng - moving.last.lng);
  moving.last = e.latlng;
  moving.moved = true;
});
let movedAt = 0;
function endMove() {
  if (!moving) return;
  const sh = moving.sh, moved = moving.moved;
  moving = null;
  map.dragging.enable();
  commitEdit();
  if (moved) { movedAt = Date.now(); refreshShape(sh); renderAll(); save(); }
  buildHandles(sh);
}
map.on('mouseup', endMove);
document.addEventListener('mouseup', endMove);

/* ------------------------------------------------------------------ *
 * 7. Drawing
 * ------------------------------------------------------------------ */
let tool = null;          // 'area' | 'line' | null
let draft = [];
let ghost = null;
let shiftHeld = false;

function setTool(t) {
  cancelDraw();
  tool = (tool === t) ? null : t;
  document.querySelectorAll('.tool').forEach(b => b.classList.toggle('on', b.dataset.tool === tool));
  map.getContainer().classList.toggle('drawing', !!tool);
  map.doubleClickZoom[tool ? 'disable' : 'enable']();
  if (tool) { state.selected = null; clearHandles(); state.shapes.forEach(s => s.layer.setStyle(styleFor(s))); }
  renderAll();
}

function drawHint() {
  const h = $('drawHint');
  if (!tool) {
    const picked = pickedShapes().filter(s => s.kind === 'area').length;
    h.textContent = !state.shapes.length
      ? ''
      : picked >= 2 ? `${picked} areas picked — Merge joins them into one shape.`
      : picked === 1 ? 'Ctrl+click another area to pick it too, then Merge.'
      : 'Click a shape to select it. Drag inside it to move it.';
    return;
  }
  const min = tool === 'line' ? 2 : 3;
  if (draft.length === 0) {
    h.textContent = tool === 'line'
      ? 'Click the start of the run.' : 'Click each corner of the roof.';
  } else {
    h.textContent = `${draft.length} point${draft.length > 1 ? 's' : ''} — ` +
      (draft.length >= min ? 'double-click or press Enter to finish. ' : '') +
      'Hold Shift for right angles.';
  }
}

function snapTo(cur) {
  if (!shiftHeld || draft.length === 0) return cur;
  const prev = draft[draft.length - 1];
  const prev2 = draft.length >= 2 ? draft[draft.length - 2] : null;
  const f = frameAt(prev[0], prev[1]);
  const c = f.toXY(cur[0], cur[1]);
  const len = Math.hypot(c[0], c[1]);
  if (len < 1e-6) return cur;
  let base = 0, step = Math.PI / 4;
  if (prev2) {
    const p2 = f.toXY(prev2[0], prev2[1]);
    base = Math.atan2(-p2[1], -p2[0]);
    step = Math.PI / 2;
  }
  const ang = base + Math.round((Math.atan2(c[1], c[0]) - base) / step) * step;
  return f.toLL(Math.cos(ang) * len, Math.sin(ang) * len);
}

function drawGhost(pts) {
  if (ghost) { map.removeLayer(ghost); ghost = null; }
  if (!pts || pts.length < 2) return;
  const opt = { color: '#fbbf24', weight: 2, dashArray: '5 4', fillColor: '#fbbf24', fillOpacity: 0.15 };
  ghost = (tool === 'area' && pts.length >= 3 ? L.polygon(pts, opt) : L.polyline(pts, opt)).addTo(map);
}

function commitDraft(pts) {
  const min = tool === 'line' ? 2 : 3;
  if (!pts || pts.length < min) return;
  const kind = tool;
  cancelDraw();
  pushUndo();
  const sh = addShape(pts, { kind });
  renderAll(); save();
  setTool(null);
  select(sh.id);
}

function cancelDraw() {
  draft = [];
  if (ghost) { map.removeLayer(ghost); ghost = null; }
  setReadout('');
  drawHint();
}

map.on('click', (e) => {
  const ll = [e.latlng.lat, e.latlng.lng];
  // a body-drag ends with a click; don't let it deselect what was just moved
  if (!tool) { if (!moving && Date.now() - movedAt > 250) select(null); return; }
  if (tool === 'area' && draft.length >= 3) {
    const first = map.latLngToContainerPoint(draft[0]);
    if (first.distanceTo(map.latLngToContainerPoint(ll)) < 12) { commitDraft(draft.slice()); return; }
  }
  draft.push(snapTo(ll));
  drawHint();
});

map.on('mousemove', (e) => {
  const ll = [e.latlng.lat, e.latlng.lng];
  if (!tool) {
    if (!moving) setReadout(`${e.latlng.lat.toFixed(6)}, ${e.latlng.lng.toFixed(6)}`);
    return;
  }
  if (draft.length === 0) return;
  const pts = draft.concat([snapTo(ll)]);
  drawGhost(pts);
  liveReadout(measureOf(pts, tool === 'line' ? 'line' : 'area'), tool);
});

map.on('dblclick', (e) => {
  if (!tool) return;
  L.DomEvent.stop(e);
  const min = tool === 'line' ? 2 : 3;
  if (draft.length >= min + 1) draft.pop();   // the dblclick's second click added a duplicate
  commitDraft(draft.slice());
});

function liveReadout(m, kind) {
  if (kind === 'line') {
    const last = m.segs.length ? m.segs[m.segs.length - 1] : 0;
    setReadout(`${fmtLen(m.len)} total${m.segs.length > 1 ? ` &nbsp;·&nbsp; ${fmtLen(last)} this leg` : ''}`, true);
    return;
  }
  const d = m.dims ? ` &nbsp;·&nbsp; ${fmtLen(m.dims.long, false)} × ${fmtLen(m.dims.short)}` : '';
  setReadout(`${fmtArea(m.area)}${d}`, true);
}
function setReadout(html, isHtml) {
  const el = $('readout');
  if (isHtml) el.innerHTML = html; else el.textContent = html;
}

/* ------------------------------------------------------------------ *
 * 8. Offset
 * ------------------------------------------------------------------ */
/* What the offset makes is the ground it adds and nothing else: the selected
   shape is left exactly where it is, and the new shape covers only the strip
   outside it. Anything else would count the same ground twice in the total. */
function doOffset() {
  const sh = selectedShape();
  if (!sh) return;
  const v = parseFloat($('offsetDist').value);
  if (!(v > 0)) { $('offsetMsg').textContent = 'Enter a distance greater than zero.'; return; }
  state.offsetDist = v;
  const metres = v / LEN_UNITS[state.lenUnit].per_m;

  const sides = pickedSides(sh);
  if (sides && !sides.size) {
    $('offsetMsg').textContent = 'Pick at least one side, or select them all.';
    return;
  }
  // All the way round, the strip is a ring: the rounded-corner buffer with the
  // shape itself cut out of it as a hole. Along a few sides it is an open band
  // with nothing to cut out -- see offsetStrip.
  const partial = sides && sides.size < sh.pts.length;
  const bands = partial ? offsetStrip(sh.pts, metres, sides) : null;
  // An area gets mitred corners, so the ring stays editable. A line keeps the
  // rounded buffer: a corridor's ends are caps, and there is no corner there
  // for two straight sides to meet at.
  const outer = partial ? null
    : sh.kind === 'area' ? offsetMitred(sh.pts, metres)
                         : offsetGeometry(sh.pts, sh.kind, metres);
  if (!bands && !outer) {
    $('offsetMsg').textContent = 'That distance collapses the shape — try a smaller one.';
    return;
  }

  const tag = `+${num(v, 0)} ${LEN_UNITS[state.lenUnit].label}`;
  const made = [];
  pushUndo();
  if (partial) {
    const list = [...sides].sort((a, b) => a - b).map(i => i + 1);
    const which = list.length <= 4
      ? `side${list.length === 1 ? '' : 's'} ${list.join(', ')}` : `${list.length} sides`;
    bands.forEach((band, i) => made.push(addShape(band, { kind: 'area',
      name: `${sh.name} ${tag} (${which}${bands.length > 1 ? `, part ${i + 1}` : ''})` })));
  } else {
    // A line encloses no ground of its own, so its offset is the whole
    // corridor and there is nothing to cut out of it.
    made.push(addShape(outer, { kind: 'area', name: `${sh.name} ${tag}`,
                                holes: sh.kind === 'area' ? [sh.pts.slice()] : null }));
  }

  renderAll(); save();
  select(made[made.length - 1].id);
  const added = made.reduce((a, m) => a + m.m.area, 0);
  offsetNote = sh.kind === 'line'
    ? `New corridor ${fmtLen(metres * 2)} wide — ${fmtArea(added)}.`
    : `New strip is ${fmtArea(added)}` +
      (made.length > 1 ? `, in ${made.length} pieces` : '') +
      `. “${escapeHtml(sh.name)}” is unchanged.`;
  renderOffsetBox();
}

/* ------------------------------------------------------------------ *
 * 9. Sidebar rendering
 * ------------------------------------------------------------------ */
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

function renderAll() { renderList(); renderTotals(); drawHint(); renderOffsetBox(); renderMerge(); }

function renderMerge() {
  const areas = pickedShapes().filter(s => s.kind === 'area');
  const b = $('mergeBtn');
  b.hidden = areas.length < 2;
  b.textContent = `Merge ${areas.length}`;
  b.title = `Join ${areas.length} areas into one shape`;
}

/* Joins the picked areas into as few shapes as they actually make. Ground
   under two of them at once is counted once afterwards, which is the whole
   point -- it is why this is a union and not just a shape with more corners. */
function doMerge() {
  const picked = pickedShapes();
  const areas = picked.filter(s => s.kind === 'area');
  const lines = picked.length - areas.length;
  if (areas.length < 2) return;

  const merged = unionShapes(areas.map(s => ({ pts: s.pts, holes: s.holes })));
  if (!merged) { ioMsg('Could not merge those. Check none of them crosses itself.'); return; }
  const before = areas.reduce((a, s) => a + measureShape(s).area, 0);
  const after = merged.reduce((a, r) => a + measureArea(r.pts, r.holes).area, 0);
  const overlap = before - after;
  // Shapes that only meet along an edge -- a ring and the building it was
  // measured around -- overlap by nothing at all, and joining them is still
  // worth doing. What says nothing happened is coming back with as many
  // shapes as went in and the same ground covered.
  if (merged.length >= areas.length && overlap < before * 1e-9) {
    ioMsg('Those areas do not touch, so joining them would not change anything.');
    return;
  }

  pushUndo();
  const first = areas[0];
  const base = `${first.name} + ${areas.length - 1} more`;
  const color = first.color, mode = first.mode;
  areas.forEach(s => removeShape(s.id));
  const made = merged.map((r, i) => addShape(r.pts, {
    kind: 'area', holes: r.holes, color, mode,
    name: base + (merged.length > 1 ? ` · part ${i + 1}` : '')
  }));
  renderAll(); save();
  select(made[0].id);
  ioMsg(`Merged ${areas.length} areas into ${made.length} — ${fmtArea(after)}.` +
        (overlap > before * 1e-9 ? ` ${fmtArea(overlap)} of overlap now counted once.` : '') +
        (lines ? ` ${lines} line${lines > 1 ? 's were' : ' was'} left alone.` : ''));
}

/* Which sides of the selected polygon the offset should push out. `sides` of
   null means all of them, which is the plain buffer. */
let sidePick = { id: null, sides: null };
let sideLayer = null;
let offsetNote = '';   // the result of the last offset, until the selection changes

function pickedSides(sh) {
  if (!sh || sh.kind !== 'area') return null;
  if (sidePick.id !== sh.id || !sidePick.sides) return null;
  return sidePick.sides;
}

function drawSideHighlight() {
  if (sideLayer) { map.removeLayer(sideLayer); sideLayer = null; }
  const sh = selectedShape();
  const sides = pickedSides(sh);
  if (!sh || !sides || sides.size === sh.pts.length) return;
  const segs = [];
  for (const i of sides) {
    segs.push([sh.pts[i], sh.pts[(i + 1) % sh.pts.length]]);
  }
  if (!segs.length) return;
  sideLayer = L.polyline(segs, {
    color: '#fff', weight: 7, opacity: 0.85, dashArray: '2 8', lineCap: 'round'
  }).addTo(map);
}

/* The chip numbers, out on the map at the middle of the side they name, so it
   is obvious which side is which before pressing Create. The Numbers button
   turns them on and they stay until it is pressed again; clicking one picks or
   drops that side, exactly as its chip does. */
let sideNums = [];
let hotSide = null;
let hotLayer = null;

function drawSideNums() {
  sideNums.forEach(m => map.removeLayer(m));
  sideNums = [];
  const sh = selectedShape();
  if (!state.sideNums || !sh || sh.kind !== 'area' || sh.hidden) return;
  const sides = pickedSides(sh);
  const n = sh.pts.length;
  // The midpoint of a side is already taken -- that is where the hollow "add a
  // corner" handle sits. The badge is nudged a fixed number of pixels clear of
  // it, outwards from the shape's centre, so both stay clickable at any zoom.
  const NUDGE = 20;
  const c = map.latLngToLayerPoint(centreOf(sh.pts));
  for (let i = 0; i < n; i++) {
    const a = sh.pts[i], b = sh.pts[(i + 1) % n];
    const mid = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
    const on = !sides || sides.has(i);
    const q = map.latLngToLayerPoint(mid);
    const len = Math.hypot(q.x - c.x, q.y - c.y) || 1;
    const ux = (q.x - c.x) / len, uy = (q.y - c.y) / len;
    const mk = L.marker(mid, {
      // Leaflet draws the icon at the point minus its anchor, so shrinking the
      // anchor along the outward direction pushes the badge that way.
      icon: L.divIcon({ className: 'sideNum' + (on ? ' on' : '') + (hotSide === i ? ' hot' : ''),
                        iconSize: [22, 22], html: String(i + 1),
                        iconAnchor: [11 - ux * NUDGE, 11 - uy * NUDGE] }),
      pane: 'handles', zIndexOffset: 1200, keyboard: false,
      title: `Side ${i + 1} — click to ${on ? 'leave it out' : 'put it back'}`
    }).addTo(map);
    // stop, or the map's own click handler deselects the shape underneath
    mk.on('click', (e) => { L.DomEvent.stop(e); toggleSide(i); });
    mk.on('mouseover', () => setHotSide(i));
    mk.on('mouseout', () => setHotSide(null));
    sideNums.push(mk);
  }
}

/* Take one side in or out of the offset. The first click on any side turns the
   implicit "all sides" into a real set, so that dropping one leaves the rest
   selected rather than starting from nothing. */
function toggleSide(i) {
  const sh = selectedShape();
  if (!sh || sh.kind !== 'area') return;
  const sides = sidePick.sides || new Set(sh.pts.map((_, k) => k));
  if (sides.has(i)) sides.delete(i); else sides.add(i);
  sidePick = { id: sh.id, sides };
  renderSides();
  drawSideHighlight();
  drawSideNums();
}

function drawHotSide() {
  if (hotLayer) { map.removeLayer(hotLayer); hotLayer = null; }
  const sh = selectedShape();
  if (hotSide === null || !sh || sh.kind !== 'area' || sh.hidden) return;
  const n = sh.pts.length;
  if (hotSide >= n) return;
  hotLayer = L.polyline([sh.pts[hotSide], sh.pts[(hotSide + 1) % n]], {
    color: '#fbbf24', weight: 9, opacity: 0.9, lineCap: 'round', interactive: false
  }).addTo(map);
}

/** Which side the pointer is over, or null once it leaves.
    This repaints the markers rather than rebuilding them: hovering a number
    would otherwise destroy the very element the pointer is pressing on, and
    the click that followed would land on nothing. */
function setHotSide(i) {
  if (hotSide === i) return;
  hotSide = i;
  drawHotSide();
  sideNums.forEach((m, k) => {
    const el = m.getElement();
    if (el) el.classList.toggle('hot', hotSide === k);
  });
}

/** The Numbers button. */
function setSideNums(on) {
  state.sideNums = !!on;
  if (!state.sideNums) hotSide = null;
  drawHotSide();
  drawSideNums();
  renderOffsetBox();
  save();
}

function renderSides() {
  const box = $('sidesBox'), chips = $('sidesChips');
  const sh = selectedShape();
  if (!sh || sh.kind !== 'area' || sh.hidden) { box.hidden = true; chips.innerHTML = ''; return; }

  box.hidden = false;
  if (sidePick.id !== sh.id) sidePick = { id: sh.id, sides: null };
  const lens = sideLengths(sh.pts);
  chips.innerHTML = '';
  lens.forEach((len, i) => {
    const on = !sidePick.sides || sidePick.sides.has(i);
    const b = document.createElement('button');
    b.className = 'chip' + (on ? ' on' : '');
    b.textContent = String(i + 1);
    b.title = `Side ${i + 1} — ${fmtLen(len)}`;
    b.addEventListener('click', () => toggleSide(i));
    b.addEventListener('mouseenter', () => { b.classList.add('hot'); setHotSide(i); });
    b.addEventListener('mouseleave', () => { b.classList.remove('hot'); setHotSide(null); });
    chips.appendChild(b);
  });

  const picked = sidePick.sides ? sidePick.sides.size : lens.length;
  $('sidesCount').textContent = picked === lens.length
    ? `all ${lens.length} sides`
    : `${picked} of ${lens.length} sides`;
}

function renderOffsetBox() {
  const sh = selectedShape();
  renderSides();
  $('offsetBtn').disabled = !sh || sh.hidden;
  $('sidesNums').classList.toggle('on', !!state.sideNums);
  $('sidesNums').setAttribute('aria-pressed', String(!!state.sideNums));
  $('offsetUnit').textContent = LEN_UNITS[state.lenUnit].label;
  $('captureUnit').textContent = LEN_UNITS[state.lenUnit].label;
  if (offsetNote) { $('offsetMsg').innerHTML = offsetNote; return; }
  if (!sh) { $('offsetMsg').textContent = 'Select a shape to offset it.'; return; }
  const sides = pickedSides(sh);
  $('offsetMsg').textContent = sides
    ? `Makes a new strip along ${sides.size} side${sides.size === 1 ? '' : 's'} of “${sh.name}”.`
    : `Makes a new strip right around “${sh.name}”, with the shape itself left out of it.`;
}

function renderList() {
  const box = $('list');
  box.innerHTML = '';
  const hiddenCount = state.shapes.filter(s => s.hidden).length;
  $('listHead').hidden = !state.shapes.length;
  $('listCount').textContent = state.shapes.length
    ? `${state.shapes.length} shape${state.shapes.length > 1 ? 's' : ''}` +
      (hiddenCount ? ` · ${hiddenCount} hidden` : '')
    : '';
  for (const sh of state.shapes) {
    const m = measureShape(sh);
    const isLine = sh.kind === 'line';
    const el = document.createElement('div');
    el.className = 'item' + (state.selected === sh.id || multi.has(sh.id) ? ' sel' : '') +
                   (sh.mode === 'subtract' && !isLine ? ' subtract' : '') + (isLine ? ' line' : '') +
                   (sh.hidden ? ' hid' : '');
    el.innerHTML =
      `<div class="top">
         <button class="eye" title="${sh.hidden ? 'Show on the map and count it' : 'Hide it and leave it out of the total'}">${sh.hidden ? '&#128584;' : '&#128065;'}</button>
         <input type="color" class="swatch" value="${colorOf(sh)}" title="Change colour">
         <input class="nm" value="${escapeHtml(sh.name)}" spellcheck="false">
         ${isLine ? '' : `<button class="pm" title="Add to / subtract from the total">${sh.mode === 'subtract' ? '&minus;' : '+'}</button>`}
         <button class="x" title="Delete">&times;</button>
       </div>
       <div class="big">${isLine ? fmtLen(m.len)
          : (sh.mode === 'subtract' ? '&minus;' : '') + fmtArea(m.area)}</div>
       <div class="meta">${isLine
          ? `${m.segs.length} segment${m.segs.length > 1 ? 's' : ''}${m.segs.length > 1 ? ' · longest ' + fmtLen(Math.max(...m.segs)) : ''}`
          : `${m.dims ? `${fmtLen(m.dims.long, false)} &times; ${fmtLen(m.dims.short)} &nbsp;·&nbsp; ` : ''}${fmtLen(m.perim)} perimeter`}</div>`;

    el.addEventListener('click', (e) => {
      if (e.target.closest('button') || e.target.closest('input')) return;
      if (e.ctrlKey || e.metaKey) toggleMulti(sh.id); else select(sh.id);
    });
    const nm = el.querySelector('.nm');
    // One undo step per visit to the field, not one per keystroke.
    nm.addEventListener('focus', beginEdit);
    nm.addEventListener('blur', commitEdit);
    nm.addEventListener('input', (e) => {
      sh.name = e.target.value; refreshShape(sh); save();
    });
    const sw = el.querySelector('.swatch');
    sw.addEventListener('focus', beginEdit);
    sw.addEventListener('change', commitEdit);
    sw.addEventListener('input', (e) => { sh.color = e.target.value; refreshShape(sh); save(); });
    sw.addEventListener('contextmenu', (e) => {      // right-click restores the default
      e.preventDefault();
      pushUndo();
      sh.color = null; refreshShape(sh); renderAll(); save();
    });
    const pm = el.querySelector('.pm');
    if (pm) pm.addEventListener('click', () => {
      pushUndo();
      sh.mode = sh.mode === 'add' ? 'subtract' : 'add';
      refreshShape(sh); renderAll(); save();
    });
    el.querySelector('.eye').addEventListener('click', () => setHidden(sh, !sh.hidden));
    el.querySelector('.x').addEventListener('click', () => { pushUndo(); removeShape(sh.id); });
    box.appendChild(el);
  }
}

function totals() {
  let add = 0, sub = 0, perim = 0, lineLen = 0, lines = 0;
  for (const sh of state.shapes) {
    if (sh.hidden) continue;
    const m = measureShape(sh);
    if (sh.kind === 'line') { lineLen += m.len; lines++; continue; }
    if (sh.mode === 'subtract') sub += m.area;
    else { add += m.area; perim += m.perim; }
  }
  return { add, sub, net: add - sub, perim, lineLen, lines };
}

function renderTotals() {
  const t = totals();
  const pf = pitchFactor();
  const roof = t.net * pf;
  const u = AREA_UNITS[state.areaUnit];
  const areas = state.shapes.filter(s => s.kind !== 'line' && !s.hidden).length;
  const hiddenCount = state.shapes.filter(s => s.hidden).length;
  let html = `<div class="lab">Total${areas ? ` &middot; ${areas} shape${areas > 1 ? 's' : ''}` : ''}</div>
    <div class="num">${fmtArea(roof)}</div>`;
  const bits = [];
  if (t.sub > 0) bits.push(`${fmtArea(t.add, false)} less ${fmtArea(t.sub)} excluded`);
  if (pf > 1.0001) bits.push(`footprint ${fmtArea(t.net)} &times; ${pf.toFixed(4)} pitch`);
  if (state.areaUnit !== 'acre' && roof > 0) bits.push(`${num(roof / 4046.8564224, 2)} acres`);
  if (t.perim > 0) bits.push(`${fmtLen(t.perim)} total perimeter`);
  if (t.lines > 0) bits.push(`${fmtLen(t.lineLen)} across ${t.lines} line${t.lines > 1 ? 's' : ''}`);
  if (hiddenCount) bits.push(`${hiddenCount} hidden, not counted`);
  if (bits.length) html += `<div class="sub">${bits.join('<br>')}</div>`;
  if (state.rate > 0 && roof > 0) {
    const cost = roof * u.per_m2 * state.rate;
    html += `<div class="sub cost">≈ ${cost.toLocaleString(undefined, { style: 'currency', currency: 'USD', maximumFractionDigits: 0 })}
             at ${state.rate}/${u.label}</div>`;
  }
  $('totals').innerHTML = html;
}

/* ------------------------------------------------------------------ *
 * 10. Export / import
 * ------------------------------------------------------------------ */
const nativeApi = () => (window.pywebview && window.pywebview.api) || null;

function download(name, mime, text) {
  const url = URL.createObjectURL(new Blob([text], { type: mime }));
  const a = document.createElement('a');
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}

async function saveText(name, mime, text, what) {
  const api = nativeApi();
  if (api) {
    const path = await api.save_file(name, text);
    ioMsg(path ? `Saved to ${path}` : 'Save cancelled.');
    return;
  }
  download(name, mime, text);
  ioMsg(what);
}

async function saveImage(name, canvas, what) {
  const api = nativeApi();
  if (api && api.save_image) {
    const b64 = canvas.toDataURL('image/png').split(',')[1];
    const path = await api.save_image(name, b64);
    ioMsg(path ? `Saved to ${path}` : 'Save cancelled.');
    return;
  }
  const blob = await new Promise(res => canvas.toBlob(res, 'image/png'));
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
  ioMsg(what);
}

const stamp = () => new Date().toISOString().slice(0, 10);

/** #rrggbb -> KML's aabbggrr byte order. */
function kmlColor(hex, aa) {
  const h = hex.replace('#', '');
  return aa + h.slice(4, 6) + h.slice(2, 4) + h.slice(0, 2);
}

function toKml() {

  const marks = state.shapes.map(sh => {
    const m = measureShape(sh);
    const isLine = sh.kind === 'line';
    const col = colorOf(sh);
    const coords = (r) => r.concat(isLine ? [] : [r[0]])
      .map(p => `${p[1].toFixed(9)},${p[0].toFixed(9)},0`).join(' ');
    const inner = (sh.holes || []).map(h =>
      `<innerBoundaryIs><LinearRing><coordinates>${coords(h)}</coordinates></LinearRing></innerBoundaryIs>`).join('\n      ');
    const desc = isLine
      ? [`Length: ${fmtLen(m.len)} (${num(m.len, 2)} m)`, `Segments: ${m.segs.length}`]
      : [`Area: ${fmtArea(m.area)} (${num(m.area, 2)} sq m)`,
         m.dims ? `Dimensions: ${fmtLen(m.dims.long, false)} x ${fmtLen(m.dims.short)}` : null,
         `Perimeter: ${fmtLen(m.perim)}`,
         sh.mode === 'subtract' ? 'Excluded from the site total.' : null];
    const geom = isLine
      ? `<LineString><tessellate>1</tessellate><coordinates>${coords(sh.pts)}</coordinates></LineString>`
      : `<Polygon><altitudeMode>clampToGround</altitudeMode>
      <outerBoundaryIs><LinearRing><coordinates>${coords(sh.pts)}</coordinates></LinearRing></outerBoundaryIs>
      ${inner}
    </Polygon>`;
    return `  <Placemark>
    <name>${escapeHtml(sh.name)}</name>
    <visibility>${sh.hidden ? 0 : 1}</visibility>
    <description><![CDATA[${desc.filter(Boolean).join('<br>')}]]></description>
    <Style>
      <LineStyle><color>${kmlColor(col, 'ff')}</color><width>2.4</width></LineStyle>
      <PolyStyle><color>${isLine ? '00000000' : kmlColor(col, '59')}</color></PolyStyle>
    </Style>
    <ExtendedData>
      <Data name="kind"><value>${sh.kind}</value></Data>
      <Data name="mode"><value>${sh.mode}</value></Data>
      <Data name="color"><value>${col}</value></Data>
      <Data name="hidden"><value>${sh.hidden ? 1 : 0}</value></Data>
      ${isLine ? `<Data name="length_m"><value>${m.len.toFixed(3)}</value></Data>`
               : `<Data name="area_sqm"><value>${m.area.toFixed(3)}</value></Data>`}
    </ExtendedData>
    ${geom}
  </Placemark>`;
  }).join('\n');

  const t = totals();
  return `<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>Site measurements ${stamp()}</name>
  <description><![CDATA[Total ${fmtArea(t.net * pitchFactor())} across ${state.shapes.length} shape(s).<br>Measured with Satellite Measurement Utility.]]></description>
${marks}
</Document>
</kml>`;
}

function toGeoJson() {
  return JSON.stringify({
    type: 'FeatureCollection',
    features: state.shapes.map(sh => {
      const m = measureShape(sh);
      const isLine = sh.kind === 'line';
      return {
        type: 'Feature',
        properties: isLine
          ? { name: sh.name, kind: 'line', hidden: !!sh.hidden,
              length_m: +m.len.toFixed(3),
              length_ft: +(m.len * 3.280839895013123).toFixed(1),
              stroke: colorOf(sh) }
          : { name: sh.name, kind: 'area', mode: sh.mode, hidden: !!sh.hidden,
              area_sqm: +m.area.toFixed(3),
              area_sqft: +(m.area * 10.763910416709722).toFixed(1),
              perimeter_m: +m.perim.toFixed(3),
              stroke: colorOf(sh), fill: colorOf(sh) },
        geometry: isLine
          ? { type: 'LineString', coordinates: sh.pts.map(p => [+p[1].toFixed(9), +p[0].toFixed(9)]) }
          : { type: 'Polygon',
              coordinates: ringsOf(sh).map(r =>
                r.concat([r[0]]).map(p => [+p[1].toFixed(9), +p[0].toFixed(9)])) }
      };
    })
  }, null, 2);
}

function parseCoordString(s) {
  return s.trim().split(/\s+/).map(tok => {
    const a = tok.split(',');
    return [parseFloat(a[1]), parseFloat(a[0])];
  }).filter(p => Number.isFinite(p[0]) && Number.isFinite(p[1]));
}

function dedupeRing(pts, closed) {
  const out = [];
  for (const p of pts) {
    const last = out[out.length - 1];
    if (!last || Math.abs(last[0] - p[0]) > 1e-12 || Math.abs(last[1] - p[1]) > 1e-12) out.push(p);
  }
  if (closed && out.length > 1) {
    const f = out[0], l = out[out.length - 1];
    if (Math.abs(f[0] - l[0]) < 1e-12 && Math.abs(f[1] - l[1]) < 1e-12) out.pop();
  }
  return out;
}

function importKml(text) {
  const doc = new DOMParser().parseFromString(text, 'application/xml');
  if (doc.getElementsByTagName('parsererror').length) throw new Error('That file is not valid KML.');
  const found = [];
  for (const pm of doc.getElementsByTagNameNS('*', 'Placemark')) {
    const nameEl = pm.getElementsByTagNameNS('*', 'name')[0];
    const base = nameEl ? nameEl.textContent.trim() : 'Imported';
    let mode = 'add';
    for (const d of pm.getElementsByTagNameNS('*', 'Data')) {
      if (d.getAttribute('name') === 'mode') {
        const v = d.getElementsByTagNameNS('*', 'value')[0];
        if (v && v.textContent.trim() === 'subtract') mode = 'subtract';
      }
    }
    // <visibility> is KML's own way of saying this, so a file from Google
    // Earth arrives with its hidden placemarks already hidden.
    let hidden = false;
    const vis = pm.getElementsByTagNameNS('*', 'visibility')[0];
    if (vis && vis.textContent.trim() === '0') hidden = true;

    const styleUrl = pm.getElementsByTagNameNS('*', 'styleUrl')[0];
    if (styleUrl && styleUrl.textContent.trim() === '#subtract') mode = 'subtract';
    let color = null;
    for (const d of pm.getElementsByTagNameNS('*', 'Data')) {
      if (d.getAttribute('name') === 'color') {
        const v = d.getElementsByTagNameNS('*', 'value')[0];
        if (v && /^#[0-9a-f]{6}$/i.test(v.textContent.trim())) color = v.textContent.trim();
      }
    }

    const polys = [];
    for (const poly of pm.getElementsByTagNameNS('*', 'Polygon')) {
      const ob = poly.getElementsByTagNameNS('*', 'outerBoundaryIs')[0] || poly;
      const co = ob.getElementsByTagNameNS('*', 'coordinates')[0];
      if (!co) continue;
      // A ring drawn by the offset comes back as a ring, not as a solid blob
      // with its hole quietly filled in.
      const holes = [];
      for (const ib of poly.getElementsByTagNameNS('*', 'innerBoundaryIs')) {
        const ic = ib.getElementsByTagNameNS('*', 'coordinates')[0];
        if (ic) holes.push(parseCoordString(ic.textContent));
      }
      polys.push({ outer: parseCoordString(co.textContent), holes });
    }
    const lines = [];
    for (const ls of pm.getElementsByTagNameNS('*', 'LineString')) {
      const co = ls.getElementsByTagNameNS('*', 'coordinates')[0];
      if (co) lines.push(parseCoordString(co.textContent));
    }
    polys.forEach((r, i) => {
      const pts = dedupeRing(r.outer, true);
      const holes = r.holes.map(h => dedupeRing(h, true)).filter(h => h.length >= 3);
      if (pts.length >= 3) found.push({ pts, holes, kind: 'area', mode, color, hidden, name: polys.length > 1 ? `${base} ${i + 1}` : base });
    });
    lines.forEach((r, i) => {
      const pts = dedupeRing(r, false);
      if (pts.length >= 2) found.push({ pts, kind: 'line', mode: 'add', color, hidden, name: lines.length > 1 ? `${base} ${i + 1}` : base });
    });
  }
  return found;
}

function importGeoJson(text) {
  const g = JSON.parse(text);
  const feats = g.type === 'FeatureCollection' ? g.features
              : g.type === 'Feature' ? [g] : [{ type: 'Feature', properties: {}, geometry: g }];
  const found = [];
  feats.forEach((f, n) => {
    const geo = f.geometry;
    if (!geo) return;
    const props = f.properties || {};
    const nm = props.name || props.Name || `Imported ${n + 1}`;
    if (geo.type === 'LineString' || geo.type === 'MultiLineString') {
      const runs = geo.type === 'LineString' ? [geo.coordinates] : geo.coordinates;
      runs.forEach((r, i) => {
        const pts = dedupeRing(r.map(c => [c[1], c[0]]), false);
        if (pts.length >= 2) found.push({ pts, kind: 'line', mode: 'add', color: props.stroke || null,
                                          hidden: !!props.hidden,
                                          name: nm + (runs.length > 1 ? ` ${i + 1}` : '') });
      });
      return;
    }
    const polys = geo.type === 'Polygon' ? [geo.coordinates]
                : geo.type === 'MultiPolygon' ? geo.coordinates : [];
    polys.forEach((poly, i) => {
      const pts = dedupeRing(poly[0].map(c => [c[1], c[0]]), true);
      const holes = poly.slice(1)
        .map(r => dedupeRing(r.map(c => [c[1], c[0]]), true)).filter(h => h.length >= 3);
      if (pts.length >= 3) {
        found.push({ pts, holes, kind: 'area', name: nm + (polys.length > 1 ? ` ${i + 1}` : ''),
                     color: props.fill || props.stroke || null,
                     hidden: !!props.hidden,
                     mode: props.mode === 'subtract' ? 'subtract' : 'add' });
      }
    });
  });
  return found;
}

function loadFound(found, label) {
  if (!found.length) { ioMsg('No shapes found in that file.'); return; }
  pushUndo();
  found.forEach(f => addShape(f.pts, { name: f.name, mode: f.mode, kind: f.kind,
                                      color: f.color, hidden: f.hidden, holes: f.holes }));
  const all = state.shapes.filter(s => !s.hidden).flatMap(s => s.pts);
  if (all.length) map.fitBounds(L.latLngBounds(all).pad(0.15));
  renderAll(); save();
  ioMsg(`Loaded ${found.length} shape${found.length > 1 ? 's' : ''} from ${label}.`);
}

function ioMsg(t) {
  $('ioMsg').textContent = t;
  clearTimeout(ioMsg._t);
  ioMsg._t = setTimeout(() => { $('ioMsg').textContent = ''; }, 6000);
}

function ingest(name, text) {
  try {
    const found = /\.kml$/i.test(name) || /<kml/i.test(text.slice(0, 400))
      ? importKml(text) : importGeoJson(text);
    loadFound(found, name);
  } catch (err) {
    ioMsg('Could not read that file: ' + err.message);
  }
}

function readFile(f) {
  const fr = new FileReader();
  fr.onload = () => ingest(f.name, String(fr.result));
  fr.readAsText(f);
}

/* ------------------------------------------------------------------ *
 * 10b. Capture
 *
 * Draws every shape that is currently shown, plus a margin of open ground
 * around them, into one PNG. The tiles are fetched again rather than scraped
 * off the map: what is on screen is whatever the window happens to be showing
 * at whatever zoom, and half a building at the edge of the viewport is not a
 * measurement anybody can use.
 *
 * Imagery for the capture is always Esri. It is the keyless layer that is
 * always available, it serves tiles with CORS so a canvas can read them back,
 * and it is licensed for this. Google's Map Tiles API forbids storing its
 * tiles, and a PNG on someone's disk is storing them.
 * ------------------------------------------------------------------ */
const CAPTURE_TILES = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}';
const CAPTURE_ATTRIB = 'Imagery © Esri, Maxar, Earthstar Geographics';
const CAPTURE_MAX_PX = 4600;   // per side; ~85 tiles at the widest
const CAPTURE_MAX_Z = 20;      // Esri's deepest imagery
const TILE = 256;

const lon2px = (lon, z) => (lon + 180) / 360 * Math.pow(2, z) * TILE;
const lat2px = (lat, z) => {
  const s = Math.sin(Math.max(-85.05, Math.min(85.05, lat)) * D2R);
  return (0.5 - Math.log((1 + s) / (1 - s)) / (4 * Math.PI)) * Math.pow(2, z) * TILE;
};

/** The area to capture: everything shown, plus `marginM` metres of ground. */
function captureBounds(marginM) {
  const shown = state.shapes.filter(s => !s.hidden && s.pts.length);
  if (!shown.length) return null;

  const pts = shown.flatMap(s => s.pts);
  let latMin = 90, latMax = -90, lonMin = 180, lonMax = -180;
  for (const p of pts) {
    latMin = Math.min(latMin, p[0]); latMax = Math.max(latMax, p[0]);
    lonMin = Math.min(lonMin, p[1]); lonMax = Math.max(lonMax, p[1]);
  }
  // Grow the box in metres, in a frame at its own centre, so the margin is the
  // same distance on every side rather than the same number of degrees.
  const f = frameAt((latMin + latMax) / 2, (lonMin + lonMax) / 2);
  const sw = f.toXY(latMin, lonMin), ne = f.toXY(latMax, lonMax);
  const a = f.toLL(sw[0] - marginM, sw[1] - marginM);
  const b = f.toLL(ne[0] + marginM, ne[1] + marginM);
  return { south: Math.min(a[0], b[0]), north: Math.max(a[0], b[0]),
           west: Math.min(a[1], b[1]), east: Math.max(a[1], b[1]), shapes: shown };
}

/** The deepest zoom whose image still fits inside CAPTURE_MAX_PX. */
function captureZoom(b) {
  for (let z = CAPTURE_MAX_Z; z > 0; z--) {
    const w = lon2px(b.east, z) - lon2px(b.west, z);
    const h = lat2px(b.south, z) - lat2px(b.north, z);
    if (w <= CAPTURE_MAX_PX && h <= CAPTURE_MAX_PX) return z;
  }
  return 1;
}

function loadTile(url) {
  return new Promise(resolve => {
    const img = new Image();
    // Without this the canvas is tainted and toDataURL throws instead of
    // producing a file.
    img.crossOrigin = 'anonymous';
    img.onload = () => resolve(img);
    img.onerror = () => resolve(null);
    img.src = url;
  });
}

function roundRect(g, x, y, w, h, r) {
  g.beginPath();
  g.moveTo(x + r, y);
  g.arcTo(x + w, y, x + w, y + h, r);
  g.arcTo(x + w, y + h, x, y + h, r);
  g.arcTo(x, y + h, x, y, r);
  g.arcTo(x, y, x + w, y, r);
  g.closePath();
}

/** Draws the capture and hands back the canvas. */
async function renderCapture(marginM, onProgress) {
  const b = captureBounds(marginM);
  if (!b) return null;

  const z = captureZoom(b);
  const x0 = lon2px(b.west, z), y0 = lat2px(b.north, z);
  const w = Math.max(1, Math.round(lon2px(b.east, z) - x0));
  const h = Math.max(1, Math.round(lat2px(b.south, z) - y0));

  // Everything drawn on top is sized off the image, not off the screen: a
  // 3,000px capture annotated at 13px is a photograph with specks on it.
  const k = Math.max(1, Math.min(3.5, Math.max(w, h) / 1300));
  const bar = Math.round(26 * k);

  const canvas = document.createElement('canvas');
  canvas.width = w; canvas.height = h + bar;
  const g = canvas.getContext('2d');
  g.fillStyle = '#0b0f14';
  g.fillRect(0, 0, canvas.width, canvas.height);

  const tx0 = Math.floor(x0 / TILE), tx1 = Math.floor((x0 + w - 1) / TILE);
  const ty0 = Math.floor(y0 / TILE), ty1 = Math.floor((y0 + h - 1) / TILE);
  const jobs = [];
  for (let ty = ty0; ty <= ty1; ty++) {
    for (let tx = tx0; tx <= tx1; tx++) {
      jobs.push({ tx, ty });
    }
  }

  let done = 0, missing = 0;
  // Eight at a time: enough to keep the connection busy, few enough that the
  // tile server does not start refusing them.
  const queue = jobs.slice();
  const workers = Array.from({ length: 8 }, async () => {
    for (let job = queue.shift(); job; job = queue.shift()) {
      const url = CAPTURE_TILES.replace('{z}', z).replace('{x}', job.tx).replace('{y}', job.ty);
      const img = await loadTile(url);
      if (img) g.drawImage(img, job.tx * TILE - x0, job.ty * TILE - y0);
      else missing++;
      if (onProgress) onProgress(++done, jobs.length);
    }
  });
  await Promise.all(workers);

  // Shapes, in the order they are listed, over the imagery.
  const px = (p) => [lon2px(p[1], z) - x0, lat2px(p[0], z) - y0];
  for (const sh of b.shapes) {
    const c = colorOf(sh);
    g.save();
    g.lineJoin = 'round';
    g.lineWidth = 3 * k;
    g.strokeStyle = c;
    if (sh.mode === 'subtract' && sh.kind !== 'line') g.setLineDash([9 * k, 6 * k]);
    g.beginPath();
    for (const r of ringsOf(sh)) {
      r.map(px).forEach((p, i) => i ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1]));
      if (sh.kind !== 'line') g.closePath();
    }
    if (sh.kind !== 'line') {
      g.fillStyle = c;
      g.globalAlpha = 0.22;
      // even-odd, so a ring's hole is left as ground rather than filled in
      // whichever way round its two boundaries happen to run.
      g.fill('evenodd');
      g.globalAlpha = 1;
    }
    g.stroke();
    g.restore();
  }

  // Labels last, so no outline is drawn over them.
  if (state.labels) {
    const fs = Math.round(13 * k);
    g.font = `600 ${fs}px system-ui, "Segoe UI", Roboto, sans-serif`;
    g.textAlign = 'center';
    g.textBaseline = 'middle';
    for (const sh of b.shapes) {
      const m = measureShape(sh);
      const has = sh.kind === 'line' ? m.len > 0 : m.area > 0;
      if (!has) continue;
      const lines = [sh.name, sh.kind === 'line' ? fmtLen(m.len)
        : (sh.mode === 'subtract' ? '−' : '') + fmtArea(m.area)];
      const at = px(m.anchor);
      const wide = Math.max(...lines.map(t => g.measureText(t).width));
      const lh = Math.round(fs * 1.3);
      const bw = wide + 16 * k, bh = lines.length * lh + 8 * k;
      g.fillStyle = 'rgba(8,12,18,0.78)';
      roundRect(g, at[0] - bw / 2, at[1] - bh / 2, bw, bh, 5 * k);
      g.fill();
      g.fillStyle = '#eef4fb';
      lines.forEach((t, i) => g.fillText(t, at[0], at[1] - bh / 2 + 4 * k + lh / 2 + i * lh));
    }
  }

  // Bottom strip: what this is, how big the margin was, and whose imagery.
  g.fillStyle = '#0b0f14';
  g.fillRect(0, h, w, bar);
  g.font = `${Math.round(12 * k)}px system-ui, "Segoe UI", Roboto, sans-serif`;
  g.textBaseline = 'middle';
  g.textAlign = 'left';
  g.fillStyle = '#9fb0c2';
  const unit = LEN_UNITS[state.lenUnit];
  g.fillText(`${new Date().toLocaleDateString()} · ${b.shapes.length} shape` +
             `${b.shapes.length === 1 ? '' : 's'} · ${num(marginM * unit.per_m, 0)} ` +
             `${unit.label} margin`, 10 * k, h + bar / 2);
  g.textAlign = 'right';
  g.fillText(CAPTURE_ATTRIB, w - 10 * k, h + bar / 2);

  // Scale bar, sized to a round number of the chosen unit.
  const mPerPx = (b.east - b.west) * 111320 * Math.cos((b.north + b.south) / 2 * D2R) / w;
  const targets = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
  const wantPx = Math.min(180 * k, w * 0.25);
  let pick = targets[0];
  for (const t of targets) if (t / unit.per_m / mPerPx <= wantPx) pick = t;
  const barPx = pick / unit.per_m / mPerPx;
  const bx = 12 * k, by = h - 20 * k;
  g.fillStyle = 'rgba(8,12,18,0.72)';
  roundRect(g, bx - 6 * k, by - 15 * k, barPx + 12 * k, 26 * k, 4 * k);
  g.fill();
  g.strokeStyle = '#eef4fb';
  g.lineWidth = 2 * k;
  g.beginPath();
  g.moveTo(bx, by); g.lineTo(bx + barPx, by);
  g.moveTo(bx, by - 4 * k); g.lineTo(bx, by + 4 * k);
  g.moveTo(bx + barPx, by - 4 * k); g.lineTo(bx + barPx, by + 4 * k);
  g.stroke();
  g.fillStyle = '#eef4fb';
  g.textAlign = 'center';
  g.font = `${Math.round(11 * k)}px system-ui, "Segoe UI", Roboto, sans-serif`;
  g.fillText(`${num(pick, 0)} ${unit.label}`, bx + barPx / 2, by - 8 * k);

  return { canvas, zoom: z, width: w, height: h, missing, shapes: b.shapes.length };
}

async function doCapture() {
  const v = parseFloat($('captureMargin').value);
  if (!(v >= 0)) { ioMsg('Enter a margin of zero or more.'); return; }
  state.captureMargin = v;
  const metres = v / LEN_UNITS[state.lenUnit].per_m;

  const shown = state.shapes.filter(s => !s.hidden);
  if (!shown.length) { ioMsg('Nothing to capture — every shape is hidden.'); return; }

  const btn = $('captureBtn');
  btn.disabled = true;
  ioMsg('Fetching imagery…');
  try {
    const out = await renderCapture(metres, (done, total) => {
      if (done % 8 === 0 || done === total) ioMsg(`Fetching imagery… ${done}/${total} tiles`);
    });
    if (!out) { ioMsg('Nothing to capture.'); return; }
    const note = out.missing ? ` (${out.missing} tile${out.missing === 1 ? '' : 's'} did not load)` : '';
    await saveImage(`site-capture-${stamp()}.png`, out.canvas,
      `Saved ${out.width}×${out.height} at zoom ${out.zoom}${note}.`);
  } catch (e) {
    ioMsg('Could not build the image: ' + e.message);
  } finally {
    btn.disabled = false;
    save();
  }
}

/* ------------------------------------------------------------------ *
 * 11. Persistence
 * ------------------------------------------------------------------ */
let saveTimer = null;
function save() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => {
    try {
      const c = map.getCenter();
      localStorage.setItem(STORE_KEY, JSON.stringify({
        shapes: state.shapes.map(plainShape),
        view: { c: [c.lat, c.lng], z: map.getZoom() },
        areaUnit: state.areaUnit, lenUnit: state.lenUnit,
        pitch: state.pitch, rate: state.rate, labels: state.labels,
        sideNums: state.sideNums,
        offsetDist: state.offsetDist, captureMargin: state.captureMargin,
        gkey: state.gkey
      }));
      $('saveState').textContent = 'Saved in this browser · ' +
        new Date().toLocaleTimeString([], { hour: 'numeric', minute: '2-digit' });
    } catch (e) { $('saveState').textContent = 'Could not autosave (storage blocked).'; }
  }, 400);
}

function restore() {
  let d;
  try {
    d = JSON.parse(localStorage.getItem(STORE_KEY) || localStorage.getItem(LEGACY_KEY) || 'null');
  } catch (e) { return false; }
  if (!d) return false;
  state.areaUnit = d.areaUnit || 'ft2';
  state.lenUnit = d.lenUnit || 'ft';
  state.pitch = d.pitch || 0;
  state.rate = d.rate || null;
  state.labels = d.labels !== false;
  state.sideNums = !!d.sideNums;
  state.offsetDist = d.offsetDist || 100;
  state.captureMargin = d.captureMargin ?? 150;
  state.gkey = d.gkey || '';
  $('areaUnit').value = state.areaUnit;
  $('lenUnit').value = state.lenUnit;
  $('pitch').value = state.pitch;
  $('rate').value = state.rate ?? '';
  $('showLabels').checked = state.labels;
  $('offsetDist').value = state.offsetDist;
  $('captureMargin').value = state.captureMargin;
  $('gkey').value = state.gkey;
  (d.shapes || []).forEach(s => addShape(s.pts, s));
  if (d.view) map.setView(d.view.c, d.view.z);
  return true;
}

/* ------------------------------------------------------------------ *
 * 11b. Undo
 *
 * Snapshot based. Every action that changes a shape records the whole set as
 * it was beforehand, and undo puts that back. Shapes are small plain objects,
 * so a copy of all of them costs nothing at this scale, and nothing can drift
 * out of step the way a per-action inverse eventually does.
 *
 * Drags and text fields call beginEdit/commitEdit instead, because a drag that
 * goes nowhere and a field that is clicked but not typed in must not fill the
 * history with steps that undo nothing.
 * ------------------------------------------------------------------ */
const UNDO_MAX = 60;
let undoStack = [];
let redoStack = [];
let pendingUndo = null;

function snapshot() {
  return JSON.stringify({
    sel: state.shapes.findIndex(s => s.id === state.selected),
    shapes: state.shapes.map(plainShape)
  });
}

function record(json) {
  if (undoStack[undoStack.length - 1] === json) return;
  undoStack.push(json);
  if (undoStack.length > UNDO_MAX) undoStack.shift();
  redoStack.length = 0;
}

/** Record the state as it is right now, before changing it. */
function pushUndo() { record(snapshot()); }

/** Remember the state before an edit that may come to nothing. */
function beginEdit() { pendingUndo = snapshot(); }

/** Keep that record, but only if the edit actually changed something. */
function commitEdit() {
  const before = pendingUndo;
  pendingUndo = null;
  if (before !== null && before !== snapshot()) record(before);
}

function applySnapshot(json) {
  const d = JSON.parse(json);
  state.shapes.forEach(sh => {
    map.removeLayer(sh.layer);
    if (sh.label) map.removeLayer(sh.label);
  });
  state.shapes = [];
  state.selected = null;
  clearHandles();
  d.shapes.forEach(sh => addShape(sh.pts, sh));
  // Ids are handed out fresh, so the selection is remembered by position.
  select(d.sel >= 0 && d.sel < state.shapes.length ? state.shapes[d.sel].id : null);
  save();
}

function undo() {
  if (!undoStack.length) { ioMsg('Nothing to undo.'); return; }
  redoStack.push(snapshot());
  applySnapshot(undoStack.pop());
  ioMsg('Undone.');
}

function redo() {
  if (!redoStack.length) { ioMsg('Nothing to redo.'); return; }
  undoStack.push(snapshot());
  applySnapshot(redoStack.pop());
  ioMsg('Redone.');
}

/** Ctrl+Z. Mid-draw the obvious thing to take back is the last point, not the
    whole shape before it -- the shape does not exist yet. */
function undoStep() {
  if (tool && draft.length) {
    draft.pop();
    drawGhost(draft);
    if (!draft.length) setReadout('');
    drawHint();
    return;
  }
  undo();
}

/* ------------------------------------------------------------------ *
 * 12. Address search
 * ------------------------------------------------------------------ */
async function doSearch() {
  const q = $('search').value.trim();
  if (!q) return;
  const msg = $('searchMsg');

  const ll = q.match(/^\s*(-?\d+(?:\.\d+)?)\s*[,\s]\s*(-?\d+(?:\.\d+)?)\s*$/);
  if (ll) {
    const lat = parseFloat(ll[1]), lng = parseFloat(ll[2]);
    if (Math.abs(lat) <= 90 && Math.abs(lng) <= 180) {
      map.setView([lat, lng], 19);
      msg.textContent = '';
      return;
    }
  }

  msg.textContent = 'Searching…';
  try {
    const url = 'https://nominatim.openstreetmap.org/search?format=jsonv2&limit=1&q=' + encodeURIComponent(q);
    const r = await fetch(url, { headers: { 'Accept': 'application/json' } });
    if (!r.ok) throw new Error('HTTP ' + r.status);
    const j = await r.json();
    if (!j.length) { msg.textContent = 'No match. Try a simpler address, or paste "lat, lng".'; return; }
    const hit = j[0];
    if (hit.boundingbox) {
      const b = hit.boundingbox.map(Number);
      map.fitBounds([[b[0], b[2]], [b[1], b[3]]], { maxZoom: 19 });
    } else {
      map.setView([+hit.lat, +hit.lon], 19);
    }
    if (map.getZoom() < 17) map.setZoom(18);
    msg.textContent = hit.display_name;
  } catch (e) {
    msg.textContent = 'Search unavailable. You can paste "lat, lng" instead.';
  }
}

/* ------------------------------------------------------------------ *
 * 13. Print sheet
 * ------------------------------------------------------------------ */
function buildPrintSheet() {
  const t = totals(), pf = pitchFactor();
  const u = AREA_UNITS[state.areaUnit];
  const rows = state.shapes.filter(sh => !sh.hidden).map(sh => {
    const m = measureShape(sh);
    if (sh.kind === 'line') {
      return `<tr><td>${escapeHtml(sh.name)} (line)</td><td>&mdash;</td>
              <td>${fmtLen(m.len)}</td><td>&mdash;</td></tr>`;
    }
    return `<tr>
      <td>${escapeHtml(sh.name)}${sh.mode === 'subtract' ? ' (excluded)' : ''}</td>
      <td>${m.dims ? `${fmtLen(m.dims.long, false)} &times; ${fmtLen(m.dims.short)}` : '&mdash;'}</td>
      <td>${fmtLen(m.perim)}</td>
      <td>${sh.mode === 'subtract' ? '&minus;' : ''}${fmtArea(m.area)}</td>
    </tr>`;
  }).join('');
  const cost = (state.rate > 0)
    ? `<tr><td colspan="3">Estimate at ${state.rate} per ${u.label}</td><td>${(t.net * pf * u.per_m2 * state.rate).toLocaleString(undefined, { style: 'currency', currency: 'USD', maximumFractionDigits: 0 })}</td></tr>`
    : '';
  $('printSheet').innerHTML = `
    <h1>Site measurements</h1>
    <div class="when">${new Date().toLocaleString()}${pf > 1.0001 ? ` &middot; roof pitch ${state.pitch}:12 applied (&times;${pf.toFixed(4)})` : ''}</div>
    <table>
      <thead><tr><th>Shape</th><th>Dimensions</th><th>Perimeter / length</th><th>Area</th></tr></thead>
      <tbody>${rows}</tbody>
      <tfoot><tr><td colspan="3">Total</td><td>${fmtArea(t.net * pf)}</td></tr>${cost}</tfoot>
    </table>
    <p class="note">Measured from satellite imagery, which is dated and may not show recent construction.
    Treat these figures as an estimate for planning, not as a survey.</p>`;
}

/* ------------------------------------------------------------------ *
 * 14. Wiring
 * ------------------------------------------------------------------ */
document.querySelectorAll('.tool').forEach(b =>
  b.addEventListener('click', () => setTool(b.dataset.tool)));
$('helpBtn').addEventListener('click', () => $('help').showModal());

$('searchBtn').addEventListener('click', doSearch);
$('search').addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); doSearch(); } });

$('offsetBtn').addEventListener('click', doOffset);
$('sidesNums').addEventListener('click', () => setSideNums(!state.sideNums));
$('sidesAll').addEventListener('click', () => {
  const sh = selectedShape();
  if (!sh) return;
  sidePick = { id: sh.id, sides: null };      // null is every side, the plain buffer
  renderAll(); drawSideHighlight();
});
$('sidesNone').addEventListener('click', () => {
  const sh = selectedShape();
  if (!sh) return;
  sidePick = { id: sh.id, sides: new Set() };
  renderAll(); drawSideHighlight();
});
$('captureBtn').addEventListener('click', doCapture);
$('captureMargin').addEventListener('input', (e) => {
  const v = parseFloat(e.target.value);
  if (v >= 0) { state.captureMargin = v; save(); }
});
$('mergeBtn').addEventListener('click', doMerge);
$('showAll').addEventListener('click', () => setAllHidden(false));
$('hideAll').addEventListener('click', () => setAllHidden(true));
$('offsetDist').addEventListener('input', (e) => {
  const v = parseFloat(e.target.value);
  if (v > 0) { state.offsetDist = v; save(); }
});

$('areaUnit').addEventListener('change', (e) => {
  state.areaUnit = e.target.value;
  state.shapes.forEach(refreshShape); renderAll(); save();
});
$('lenUnit').addEventListener('change', (e) => {
  state.lenUnit = e.target.value;
  state.shapes.forEach(refreshShape); renderAll(); save();
});
$('pitch').addEventListener('input', (e) => {
  state.pitch = Math.max(0, parseFloat(e.target.value) || 0); renderTotals(); save();
});
$('rate').addEventListener('input', (e) => {
  const v = parseFloat(e.target.value);
  state.rate = Number.isFinite(v) && v > 0 ? v : null;
  renderTotals(); save();
});
$('showLabels').addEventListener('change', (e) => {
  state.labels = e.target.checked;
  state.shapes.forEach(refreshShape); save();
});

$('gkeyApply').addEventListener('click', async () => {
  const key = $('gkey').value.trim();
  if (!key) { $('gkeyMsg').textContent = 'Paste a Google Maps Platform API key first.'; return; }
  if (await enableGoogle(key)) { state.gkey = key; save(); }
});
$('gkeyClear').addEventListener('click', () => {
  $('gkey').value = '';
  state.gkey = '';
  disableGoogle();
  if (!map.hasLayer(esri)) map.addLayer(esri);
  $('gkeyMsg').textContent = 'Back on Esri imagery.';
  save();
});

$('expKml').addEventListener('click', () => {
  if (!state.shapes.length) { ioMsg('Nothing to save yet.'); return; }
  saveText(`site-measurements-${stamp()}.kml`, 'application/vnd.google-earth.kml+xml',
           toKml(), 'Saved .kml — opens in Google Earth.');
});
$('expGeo').addEventListener('click', () => {
  if (!state.shapes.length) { ioMsg('Nothing to save yet.'); return; }
  saveText(`site-measurements-${stamp()}.geojson`, 'application/geo+json',
           toGeoJson(), 'Saved .geojson.');
});
$('impBtn').addEventListener('click', async () => {
  const api = nativeApi();
  if (api) {
    const res = await api.open_file();
    if (res) ingest(res.name, res.text);
    return;
  }
  $('impFile').click();
});
$('impFile').addEventListener('change', (e) => {
  const f = e.target.files[0];
  if (f) readFile(f);
  e.target.value = '';
});

$('printBtn').addEventListener('click', () => {
  if (!state.shapes.length) { ioMsg('Nothing to print yet.'); return; }
  buildPrintSheet();
  window.print();
});

$('clearBtn').addEventListener('click', () => {
  if (!state.shapes.length) return;
  if (!confirm(`Delete all ${state.shapes.length} shape(s)? Ctrl+Z brings them back, and a .kml keeps them for good.`)) return;
  pushUndo();
  state.shapes.slice().forEach(s => removeShape(s.id));
  renderAll(); save();
});

['dragenter', 'dragover'].forEach(ev =>
  document.addEventListener(ev, (e) => { e.preventDefault(); e.dataTransfer.dropEffect = 'copy'; }));
document.addEventListener('drop', (e) => {
  e.preventDefault();
  const f = e.dataTransfer.files[0];
  if (f) readFile(f);
});

document.addEventListener('keydown', (e) => {
  if (e.key === 'Shift') shiftHeld = true;
  const typing = /^(INPUT|SELECT|TEXTAREA)$/.test(document.activeElement.tagName);
  // While a field has focus, Ctrl+Z is the browser's own text undo, which is
  // what someone half way through renaming a shape means by it.
  if (typing || $('help').open) return;
  if (e.ctrlKey || e.metaKey) {
    const k = e.key.toLowerCase();
    if (k === 'z' && !e.shiftKey) { e.preventDefault(); undoStep(); }
    else if (k === 'y' || (k === 'z' && e.shiftKey)) { e.preventDefault(); redo(); }
    return;   // and no bare-letter shortcut fires off a browser combo
  }
  const min = tool === 'line' ? 2 : 3;
  if (e.key === 'Escape') { tool ? setTool(null) : select(null); }
  else if (e.key === 'Enter' && tool && draft.length >= min) commitDraft(draft.slice());
  else if (e.key === 'a' || e.key === 'A') setTool('area');
  else if (e.key === 'l' || e.key === 'L') setTool('line');
  else if ((e.key === 'Delete' || e.key === 'Backspace') && state.selected != null) {
    pushUndo(); removeShape(state.selected);
  }
});
document.addEventListener('keyup', (e) => { if (e.key === 'Shift') shiftHeld = false; });

map.on('moveend zoomend', save);
map.on('zoomend', drawSideNums);   // the badges are nudged in pixels, not metres
window.addEventListener('beforeprint', () => {
  if (!$('printSheet').innerHTML) buildPrintSheet();
  map.invalidateSize({ animate: false });
});
window.addEventListener('afterprint', () => map.invalidateSize({ animate: false }));

/* ------------------------------------------------------------------ *
 * 15. Boot
 * ------------------------------------------------------------------ */
if (!restore()) {
  map.setView([39.5, -98.35], 5);
  $('searchMsg').textContent = 'Start by searching for the site address.';
}
if (state.gkey) {
  $('imageryBox').open = false;
  enableGoogle(state.gkey);
}
renderAll();
