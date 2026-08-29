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

/** Area shape: area (m^2), perimeter (m), bounding dimensions, label anchor. */
function measureArea(pts) {
  const n = pts.length;
  if (n < 2) return { area: 0, perim: 0, dims: null, anchor: pts[0] || [0, 0] };
  const c = centreOf(pts);
  const f = frameAt(c[0], c[1]);
  const xy = pts.map(p => f.toXY(p[0], p[1]));
  let twice = 0, perim = 0, cx = 0, cy = 0;
  for (let i = 0; i < n; i++) {
    const a = xy[i], b = xy[(i + 1) % n];
    const cr = a[0] * b[1] - b[0] * a[1];
    twice += cr;
    cx += (a[0] + b[0]) * cr;
    cy += (a[1] + b[1]) * cr;
    perim += Math.hypot(b[0] - a[0], b[1] - a[1]);
  }
  const area = Math.abs(twice) / 2;
  const anchor = Math.abs(twice) > 1e-9 ? f.toLL(cx / (3 * twice), cy / (3 * twice)) : c;
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

const measureOf = (pts, kind) => kind === 'line' ? measureLine(pts) : measureArea(pts);

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
  const keep = cand.filter(p => distToPath(p, segs) >= d - tol);

  const outXY = [];
  for (const p of keep) {
    const last = outXY[outXY.length - 1];
    if (!last || Math.hypot(p[0] - last[0], p[1] - last[1]) > tol) outXY.push(p);
  }
  if (outXY.length > 2) {
    const a = outXY[0], b = outXY[outXY.length - 1];
    if (Math.hypot(a[0] - b[0], a[1] - b[1]) <= tol) outXY.pop();
  }
  if (outXY.length < 3) return null;
  return outXY.map(p => f.toLL(p[0], p[1]));
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
  offsetDist: 100,
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
  const sh = {
    id: nextId++,
    name: opts.name || defaultName(kind),
    kind,
    mode: opts.mode === 'subtract' ? 'subtract' : 'add',
    color: opts.color || null,
    pts,
    layer: null,
    label: null
  };
  sh.layer = (kind === 'line' ? L.polyline(pts, styleFor(sh)) : L.polygon(pts, styleFor(sh))).addTo(map);
  sh.layer.on('click', (e) => { L.DomEvent.stop(e); if (!tool) select(sh.id); });
  sh.layer.on('mousedown', (e) => startMove(sh, e));
  state.shapes.push(sh);
  refreshShape(sh);
  return sh;
}

/** A shape's own colour if it has been customised, otherwise the default
    for its kind and mode. */
const colorOf = (sh) => sh.color || (sh.kind === 'line' ? COLOR.line : COLOR[sh.mode]);

function styleFor(sh) {
  const c = colorOf(sh);
  const on = state.selected === sh.id;
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
  sh.layer.setLatLngs(sh.pts);
  sh.layer.setStyle(styleFor(sh));
  const m = measureOf(sh.pts, sh.kind);
  sh.m = m;
  if (sh.label) { map.removeLayer(sh.label); sh.label = null; }
  const has = sh.kind === 'line' ? m.len > 0 : m.area > 0;
  if (state.labels && has) {
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
  if (state.selected === id) select(null);
  else { renderAll(); save(); }
}

function select(id) {
  state.selected = id;
  clearHandles();
  state.shapes.forEach(s => s.layer.setStyle(styleFor(s)));
  const sh = state.shapes.find(s => s.id === id);
  if (sh) buildHandles(sh);
  renderAll();
}

const selectedShape = () => state.shapes.find(s => s.id === state.selected) || null;

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
    last = e.target.getLatLng();
    handles.forEach(h => { if (h !== mv) map.removeLayer(h); });
    handles = [mv];
  });
  mv.on('drag', (e) => {
    const cur = e.target.getLatLng();
    const dLat = cur.lat - last.lat, dLng = cur.lng - last.lng;
    last = cur;
    sh.pts = sh.pts.map(p => [p[0] + dLat, p[1] + dLng]);
    sh.layer.setLatLngs(sh.pts);
  });
  mv.on('dragend', () => { refreshShape(sh); renderAll(); save(); buildHandles(sh); });
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
    mk.on('drag', (e) => {
      const ll = e.target.getLatLng();
      sh.pts[i] = [ll.lat, ll.lng];
      sh.layer.setLatLngs(sh.pts);
      liveReadout(measureOf(sh.pts, sh.kind), sh.kind);
    });
    mk.on('dragend', () => { refreshShape(sh); buildHandles(sh); renderAll(); save(); setReadout(''); });
    mk.on('contextmenu', (e) => {
      L.DomEvent.stop(e);
      if (sh.pts.length <= (closed ? 3 : 2)) return;
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
  moving = { sh, last: e.latlng, moved: false };
}
map.on('mousemove', (e) => {
  if (!moving) return;
  const dLat = e.latlng.lat - moving.last.lat, dLng = e.latlng.lng - moving.last.lng;
  moving.last = e.latlng;
  moving.moved = true;
  moving.sh.pts = moving.sh.pts.map(p => [p[0] + dLat, p[1] + dLng]);
  moving.sh.layer.setLatLngs(moving.sh.pts);
});
let movedAt = 0;
function endMove() {
  if (!moving) return;
  const sh = moving.sh, moved = moving.moved;
  moving = null;
  map.dragging.enable();
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
    h.textContent = state.shapes.length
      ? 'Click a shape to select it. Drag inside it to move it.' : '';
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
function doOffset() {
  const sh = selectedShape();
  if (!sh) return;
  const v = parseFloat($('offsetDist').value);
  if (!(v > 0)) { $('offsetMsg').textContent = 'Enter a distance greater than zero.'; return; }
  state.offsetDist = v;
  const metres = v / LEN_UNITS[state.lenUnit].per_m;
  const pts = offsetGeometry(sh.pts, sh.kind, metres);
  if (!pts) {
    $('offsetMsg').textContent = 'That distance collapses the shape — try a smaller one.';
    return;
  }
  const made = addShape(pts, { name: `${sh.name} +${num(v, 0)} ${LEN_UNITS[state.lenUnit].label}`, kind: 'area' });
  renderAll(); save();
  select(made.id);
  const ring = made.m.area - (sh.kind === 'area' ? sh.m.area : 0);
  $('offsetMsg').innerHTML = sh.kind === 'area'
    ? `Ring alone is ${fmtArea(ring)}.`
    : `Corridor ${fmtLen(metres * 2)} wide.`;
}

/* ------------------------------------------------------------------ *
 * 9. Sidebar rendering
 * ------------------------------------------------------------------ */
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

function renderAll() { renderList(); renderTotals(); drawHint(); renderOffsetBox(); }

function renderOffsetBox() {
  const sh = selectedShape();
  $('offsetBtn').disabled = !sh;
  $('offsetUnit').textContent = LEN_UNITS[state.lenUnit].label;
  if (!sh) $('offsetMsg').textContent = 'Select a shape to offset it.';
  else if (!$('offsetMsg').textContent.startsWith('Ring') &&
           !$('offsetMsg').textContent.startsWith('Corridor')) {
    $('offsetMsg').textContent = `Offsets “${sh.name}” outward all the way around.`;
  }
}

function renderList() {
  const box = $('list');
  box.innerHTML = '';
  for (const sh of state.shapes) {
    const m = sh.m || measureOf(sh.pts, sh.kind);
    const isLine = sh.kind === 'line';
    const el = document.createElement('div');
    el.className = 'item' + (state.selected === sh.id ? ' sel' : '') +
                   (sh.mode === 'subtract' && !isLine ? ' subtract' : '') + (isLine ? ' line' : '');
    el.innerHTML =
      `<div class="top">
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
      select(sh.id);
    });
    el.querySelector('.nm').addEventListener('input', (e) => {
      sh.name = e.target.value; refreshShape(sh); save();
    });
    const sw = el.querySelector('.swatch');
    sw.addEventListener('input', (e) => { sh.color = e.target.value; refreshShape(sh); save(); });
    sw.addEventListener('contextmenu', (e) => {      // right-click restores the default
      e.preventDefault();
      sh.color = null; refreshShape(sh); renderAll(); save();
    });
    const pm = el.querySelector('.pm');
    if (pm) pm.addEventListener('click', () => {
      sh.mode = sh.mode === 'add' ? 'subtract' : 'add';
      refreshShape(sh); renderAll(); save();
    });
    el.querySelector('.x').addEventListener('click', () => removeShape(sh.id));
    box.appendChild(el);
  }
}

function totals() {
  let add = 0, sub = 0, perim = 0, lineLen = 0, lines = 0;
  for (const sh of state.shapes) {
    const m = sh.m || measureOf(sh.pts, sh.kind);
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
  const areas = state.shapes.filter(s => s.kind !== 'line').length;
  let html = `<div class="lab">Total${areas ? ` &middot; ${areas} shape${areas > 1 ? 's' : ''}` : ''}</div>
    <div class="num">${fmtArea(roof)}</div>`;
  const bits = [];
  if (t.sub > 0) bits.push(`${fmtArea(t.add, false)} less ${fmtArea(t.sub)} excluded`);
  if (pf > 1.0001) bits.push(`footprint ${fmtArea(t.net)} &times; ${pf.toFixed(4)} pitch`);
  if (state.areaUnit !== 'acre' && roof > 0) bits.push(`${num(roof / 4046.8564224, 2)} acres`);
  if (t.perim > 0) bits.push(`${fmtLen(t.perim)} total perimeter`);
  if (t.lines > 0) bits.push(`${fmtLen(t.lineLen)} across ${t.lines} line${t.lines > 1 ? 's' : ''}`);
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

const stamp = () => new Date().toISOString().slice(0, 10);

/** #rrggbb -> KML's aabbggrr byte order. */
function kmlColor(hex, aa) {
  const h = hex.replace('#', '');
  return aa + h.slice(4, 6) + h.slice(2, 4) + h.slice(0, 2);
}

function toKml() {

  const marks = state.shapes.map(sh => {
    const m = sh.m || measureOf(sh.pts, sh.kind);
    const isLine = sh.kind === 'line';
    const col = colorOf(sh);
    const coords = (isLine ? sh.pts : sh.pts.concat([sh.pts[0]]))
      .map(p => `${p[1].toFixed(9)},${p[0].toFixed(9)},0`).join(' ');
    const desc = isLine
      ? [`Length: ${fmtLen(m.len)} (${num(m.len, 2)} m)`, `Segments: ${m.segs.length}`]
      : [`Area: ${fmtArea(m.area)} (${num(m.area, 2)} sq m)`,
         m.dims ? `Dimensions: ${fmtLen(m.dims.long, false)} x ${fmtLen(m.dims.short)}` : null,
         `Perimeter: ${fmtLen(m.perim)}`,
         sh.mode === 'subtract' ? 'Excluded from the site total.' : null];
    const geom = isLine
      ? `<LineString><tessellate>1</tessellate><coordinates>${coords}</coordinates></LineString>`
      : `<Polygon><altitudeMode>clampToGround</altitudeMode>
      <outerBoundaryIs><LinearRing><coordinates>${coords}</coordinates></LinearRing></outerBoundaryIs>
    </Polygon>`;
    return `  <Placemark>
    <name>${escapeHtml(sh.name)}</name>
    <description><![CDATA[${desc.filter(Boolean).join('<br>')}]]></description>
    <Style>
      <LineStyle><color>${kmlColor(col, 'ff')}</color><width>2.4</width></LineStyle>
      <PolyStyle><color>${isLine ? '00000000' : kmlColor(col, '59')}</color></PolyStyle>
    </Style>
    <ExtendedData>
      <Data name="kind"><value>${sh.kind}</value></Data>
      <Data name="mode"><value>${sh.mode}</value></Data>
      <Data name="color"><value>${col}</value></Data>
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
      const m = sh.m || measureOf(sh.pts, sh.kind);
      const isLine = sh.kind === 'line';
      return {
        type: 'Feature',
        properties: isLine
          ? { name: sh.name, kind: 'line', length_m: +m.len.toFixed(3),
              length_ft: +(m.len * 3.280839895013123).toFixed(1),
              stroke: colorOf(sh) }
          : { name: sh.name, kind: 'area', mode: sh.mode,
              area_sqm: +m.area.toFixed(3),
              area_sqft: +(m.area * 10.763910416709722).toFixed(1),
              perimeter_m: +m.perim.toFixed(3),
              stroke: colorOf(sh), fill: colorOf(sh) },
        geometry: isLine
          ? { type: 'LineString', coordinates: sh.pts.map(p => [+p[1].toFixed(9), +p[0].toFixed(9)]) }
          : { type: 'Polygon',
              coordinates: [sh.pts.concat([sh.pts[0]]).map(p => [+p[1].toFixed(9), +p[0].toFixed(9)])] }
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
      if (co) polys.push(parseCoordString(co.textContent));
    }
    const lines = [];
    for (const ls of pm.getElementsByTagNameNS('*', 'LineString')) {
      const co = ls.getElementsByTagNameNS('*', 'coordinates')[0];
      if (co) lines.push(parseCoordString(co.textContent));
    }
    polys.forEach((r, i) => {
      const pts = dedupeRing(r, true);
      if (pts.length >= 3) found.push({ pts, kind: 'area', mode, color, name: polys.length > 1 ? `${base} ${i + 1}` : base });
    });
    lines.forEach((r, i) => {
      const pts = dedupeRing(r, false);
      if (pts.length >= 2) found.push({ pts, kind: 'line', mode: 'add', color, name: lines.length > 1 ? `${base} ${i + 1}` : base });
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
                                          name: nm + (runs.length > 1 ? ` ${i + 1}` : '') });
      });
      return;
    }
    const polys = geo.type === 'Polygon' ? [geo.coordinates]
                : geo.type === 'MultiPolygon' ? geo.coordinates : [];
    polys.forEach((poly, i) => {
      const pts = dedupeRing(poly[0].map(c => [c[1], c[0]]), true);
      if (pts.length >= 3) {
        found.push({ pts, kind: 'area', name: nm + (polys.length > 1 ? ` ${i + 1}` : ''),
                     color: props.fill || props.stroke || null,
                     mode: props.mode === 'subtract' ? 'subtract' : 'add' });
      }
    });
  });
  return found;
}

function loadFound(found, label) {
  if (!found.length) { ioMsg('No shapes found in that file.'); return; }
  found.forEach(f => addShape(f.pts, { name: f.name, mode: f.mode, kind: f.kind, color: f.color }));
  const all = state.shapes.flatMap(s => s.pts);
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
 * 11. Persistence
 * ------------------------------------------------------------------ */
let saveTimer = null;
function save() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => {
    try {
      const c = map.getCenter();
      localStorage.setItem(STORE_KEY, JSON.stringify({
        shapes: state.shapes.map(s => ({ name: s.name, kind: s.kind, mode: s.mode,
                                         color: s.color, pts: s.pts })),
        view: { c: [c.lat, c.lng], z: map.getZoom() },
        areaUnit: state.areaUnit, lenUnit: state.lenUnit,
        pitch: state.pitch, rate: state.rate, labels: state.labels,
        offsetDist: state.offsetDist, gkey: state.gkey
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
  state.offsetDist = d.offsetDist || 100;
  state.gkey = d.gkey || '';
  $('areaUnit').value = state.areaUnit;
  $('lenUnit').value = state.lenUnit;
  $('pitch').value = state.pitch;
  $('rate').value = state.rate ?? '';
  $('showLabels').checked = state.labels;
  $('offsetDist').value = state.offsetDist;
  $('gkey').value = state.gkey;
  (d.shapes || []).forEach(s => addShape(s.pts, { name: s.name, mode: s.mode, kind: s.kind,
                                                  color: s.color }));
  if (d.view) map.setView(d.view.c, d.view.z);
  return true;
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
  const rows = state.shapes.map(sh => {
    const m = sh.m || measureOf(sh.pts, sh.kind);
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
  if (!confirm(`Delete all ${state.shapes.length} shape(s)? Save a .kml first if you want to keep them.`)) return;
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
  if (typing || $('help').open) return;
  const min = tool === 'line' ? 2 : 3;
  if (e.key === 'Escape') { tool ? setTool(null) : select(null); }
  else if (e.key === 'Enter' && tool && draft.length >= min) commitDraft(draft.slice());
  else if (e.key === 'a' || e.key === 'A') setTool('area');
  else if (e.key === 'l' || e.key === 'L') setTool('line');
  else if ((e.key === 'Delete' || e.key === 'Backspace') && state.selected != null) removeShape(state.selected);
});
document.addEventListener('keyup', (e) => { if (e.key === 'Shift') shiftHeld = false; });

map.on('moveend zoomend', save);
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
