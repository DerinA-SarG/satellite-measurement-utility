# CLAUDE.md — Satellite Measurement Utility

Instructions for Claude working in this repo. Read before changing anything.

## What this is

A zero-dependency web app that measures areas and distances on satellite
imagery, plus a thin Python shell that wraps it as a desktop `.exe`. No
account, no API key, no backend, no build step for the web version. Everything
the browser needs is `index.html`, `app.css`, `app.js` and `vendor/`.

## Ground rules

These are the things that make the tool correct. Do not trade them away for
convenience.

1. **Never replace the geodesy with a spherical formula.** Areas are computed
   in a local East/North metre frame built from the true WGS84 radii of
   curvature at each shape's own latitude (`frameAt`), then measured with the
   shoelace formula. Turf's `area` and Google's `computeArea` both use a sphere
   of radius 6378137 and run **0.67% high at the equator, 0.12% high at 40°N,
   0.41% low at 64°N**. Swapping one in would silently make every number worse.
2. **Never fetch Google imagery from anything but the Map Tiles API.** The
   `mt0.google.com/vt` tile endpoint that circulates online violates Google's
   terms. The only sanctioned path is `createSession` then
   `tile.googleapis.com/v1/2dtiles`, which is what `enableGoogle()` does. The
   PNG capture never uses Google at all, whatever the map is showing: Google
   forbids storing its tiles, and writing one to someone's disk is storing it.
   Esri is licensed for this and serves tiles with CORS, which the capture also
   needs -- without `crossOrigin` the canvas is tainted and `toDataURL` throws
   instead of producing a file.
3. **Keep Esri as the always-available fallback.** Google is opt-in and billed
   per tile; Esri is free and keyless. Any Google failure must fall back rather
   than leave a blank map.
4. **Keep Leaflet vendored.** No CDN. `vendor/` ships with the repo so the tool
   works if unpkg is down, and `vendor/LEAFLET-LICENSE.txt` must stay — we
   redistribute `leaflet.js` under BSD-2-Clause.
5. **No secrets, no binaries, no user data in git.** The Google key lives in
   the user's `localStorage` only. `dist/` and `icon.ico` are generated and
   gitignored; the exe ships via GitHub Releases.
6. **Preserve saved work across renames.** `restore()` reads `STORE_KEY` and
   falls back to `LEGACY_KEY`. If you change the key again, add another
   fallback rather than orphaning people's shapes.

## Repo map

```
index.html      markup; the sidebar is the whole UI
app.css         styling, including the @media print stylesheet
app.js          everything, in 15 numbered sections
vendor/         Leaflet 1.9.4 + its licence
desktop.py      native-window shell: loopback server + OS file dialogs
                (save_file for text, save_image for the capture's base64 PNG)
build_exe.py    generates icon.ico, then PyInstaller --onefile --windowed
start.bat/.sh   run it in a browser
```

`app.js` sections: 1 geodesy · 1b offset (buffer, mitred ring, per-side,
strip) · 1c union · 2 units · 3 state · 4 map/imagery · 5 shapes · 6 handles ·
7 drawing · 8 offset UI · 9 sidebar · 10 export/import · 10b capture ·
11 persistence · 11b undo · 12 search · 13 print · 14 wiring · 15 boot.

Sections 1, 1b and 1c are pure maths with no DOM or Leaflet dependency — keep
them that way so they stay testable in isolation.

The offset a user gets is **mitred**, not the rounded buffer. The rounded one
is the more faithful answer and is still there — `offsetGeometry`, used for
line corridors and still validated against `A + Pd + πd²` — but it spends
seventy vertices on the corners of a warehouse, and a shape with seventy
vertices cannot be dragged into place. `offsetMitred` reuses `offsetGeometry`'s
own filter to drop corners that fold back inside a notch narrower than twice
the offset. Do not "fix" the ring back to arcs.

## Commands

```bash
python -m http.server 8123        # run the web version
pip install -r requirements.txt   # only needed to build the exe
python build_exe.py               # -> dist/SatelliteMeasurementUtility.exe
```

There is no test runner. Verification is done by executing probes in the live
page, described below. Do that rather than assuming a change is correct.

## How to verify a change

**Any change to section 1 or 1b must be re-validated numerically.** The method
is to check against a closed form, not against the previous output.

*Area* — the exact area of a lat/lon cell on the WGS84 ellipsoid has a closed
form. Compare `measureArea` against it across latitudes; error must stay at
0.0000%:

```python
Z = lambda phi: A*A*(1-E2)/2 * (math.sin(phi)/(1-E2*math.sin(phi)**2)
                + (1/(2*E))*math.log((1+E*math.sin(phi))/(1-E*math.sin(phi))))
exact = (lon2-lon1)*D2R * (Z(lat2*D2R) - Z(lat1*D2R))
```

Sanity-check the integral itself first: `2π(Z(π/2) − Z(−π/2))` must give
510,065,622 km².

*Offset* — a buffer of distance `d` around a shape has area `A + Pd + πd²`.
Check that, and separately assert that **every** output vertex sits exactly `d`
from the source via `distToPath`. Expect ≤0.02% on area (the residual is the
polygon approximation of the corner arcs) and exact on distance.

*Partial offset* — `offsetSides(pts, d, sides)` mitres its corners instead of
rounding them, which makes a rectangle an exact test. Pushing one side of a
`w × h` rectangle out by `d` gives exactly `A + wd`; two adjacent sides give
`A + (w+h)d + d²`; all four give `A + Pd + 4d²`. Expect 0.000% on each --
there are no arcs to approximate. Check a reversed (clockwise) ring too: side
`i` must always be the edge from `pts[i]` to `pts[i+1]` whichever way the ring
was drawn, because that is what the numbered chips in the sidebar point at.

*Strip* — `offsetStrip(pts, d, sides)` is what the tool actually creates: the
same geometry with the original taken out of it, so the areas above lose their
`A`. One side of a rectangle gives `wd`, two adjacent give `(w+h)d + d²`, three
give `(2w+h)d + 2d²` — one `d²` per corner where two chosen sides meet. Expect
0.000%. Walk **all** 2ⁿ−1 subsets on a rectangle, both windings: the bug this
catches is a run of sides that wraps past vertex 0, which is indexed by run
length rather than by end vertex for exactly that reason. Sides on opposite
walls must come back as two rings, not one. All sides returns null; the ring is
`offsetGeometry` plus the shape as a hole, and its area must equal the buffer
minus the original to 0.00000%, with `measureArea(buf, [pts])` matching
`Pd + πd²` inside the same 0.02%. The ring's label anchor must land in the band
— inside the outline, outside the hole — because a ring's centroid does not.

*Mitred ring* — `offsetMitred(pts, d)` is exact, so test it against the closed
form directly: for a polygon with only right angles the offset area is
`A + Pd + 4d²`, convex corners contributing `+d²` and reflex ones `−d²`. Build
the test shape **centred on the frame origin**, or a 1e-4% round-trip error
from the two different local frames swamps the result. Check an L, a wide notch,
and a notch narrower than `2d`: the narrow one must lose corners rather than
grow a loop. Assert zero self-intersections and that every output vertex is at
least `d` from the source, both times.

*Union* — `unionShapes([{pts, holes}, …])` on rectangles in one fixed frame is
exact. Two overlapping give `A₁ + A₂ − overlap`; two sharing a wall give the
sum and must come back as **four** corners, not six, because `dropCollinear`
takes the old wall ends out; one inside another gives the outer alone; disjoint
ones come back as two shapes. A U plus a cap across its mouth must return one
outline **with a hole**, not a filled rectangle — that is the case that proves
the loop tracing and hole assignment. A ring merged with the building it was
measured around must close its hole and return the plain outer rectangle.
Beware the disjointness test in `doMerge`: shapes meeting along an edge overlap
by *nothing*, so zero overlap alone must not be read as "nothing happened".

*Capture* — `captureBounds(marginM)` must sit exactly `marginM` from the
outermost shape on all four sides (check in a local frame, not in degrees), and
must ignore hidden shapes entirely. `renderCapture` reports `missing`; a
non-zero count means tiles failed and the image has holes in it.

*In the browser* — open the page and run probes with the JS tool:

```js
measureArea(pts)            // area m², perimeter, minRect dims, label anchor
offsetGeometry(pts, kind, d)
importKml(toKml())          // round-trip must preserve pts, kind, mode, colour
```

Check `read_console_messages` for errors after every change.

## Conventions

- Plain ES2020, no framework, no build. Strict mode, `const`/`let`.
- Two-space indent in JS/CSS/HTML, four in Python.
- Shapes are `{id, name, kind:'area'|'line', mode:'add'|'subtract', color,
  hidden, pts, holes}` where `pts` is `[[lat,lng],…]`, unclosed. Closing points
  are added only at export time.
- `holes` is null on nearly everything. It is a list of rings cut out of `pts`
  — how a full-ring offset holds the shape it was measured around. Rings come
  off the area and add their edge to the perimeter, Leaflet takes them as extra
  entries in `setLatLngs`, KML as `innerBoundaryIs`, GeoJSON as the rings after
  the first, and the capture fills even-odd so the hole stays as ground. Read
  a shape through `ringsOf`/`latLngsOf`/`measureShape` rather than `sh.pts`, so
  none of that gets dropped. Vertex handles stay on the outline only; a hole is
  geometry the offset put there, not something to drag.
- `hidden` is a view decision: the shape leaves the map and the totals but
  stays in the list, in the save file, and in both export formats -- KML says
  it with its own `<visibility>` tag, GeoJSON with a `hidden` property. Anything
  that sums shapes has to skip it; anything that lists them must not.
- Colour is `null` when the shape uses its kind/mode default; `colorOf()`
  resolves it. Never write the default into `color`.
- The `.exe` keeps nothing between runs, and that is deliberate: pywebview
  defaults to `private_mode`, and the loopback server takes a fresh port each
  launch so the origin changes anyway. Export is the persistence story. Do not
  "fix" it without being asked.
- Anything that changes a shape has to record an undo step first. Discrete
  actions call `pushUndo()` before the change; drags and text fields call
  `beginEdit()` when they start and `commitEdit()` when they end, which throws
  the record away if nothing actually moved. A new action that skips this does
  not break — it just quietly makes Ctrl+Z jump over it, which is worse.
- Comments explain *why*, not *what*. Match the existing density — sparse, with
  a short block above anything non-obvious.
- British spelling in UI text and comments ("colour", "centre").

## Things that will bite you

Each of these cost real debugging time. They are fixed; do not undo them.

- **pywebview introspection.** Storing the window on the `Api` object sends
  pywebview recursing through WebView2 COM objects and spraying
  `maximum recursion depth exceeded` at startup. `Api` must hold only plain
  methods; get the window via `webview.windows[0]` at call time.
- **Leaflet pane order.** Tooltips (650) sit above markers (600), so edit
  handles vanish under measurement labels. Handles live in a custom `handles`
  pane at 680.
- **Thin polylines are unhittable.** Dragging a 4px line body grabs the map and
  pans it instead. That is why every selected shape gets a centre move grip.
- **Double-click adds a duplicate vertex.** Leaflet fires `click`, `click`,
  `dblclick`. The `dblclick` handler pops the extra point.
- **Leaflet positions divIcons with `transform`.** Anything in `.vtx` or
  `.sideNum` that sets a transform of its own — a `scale()` to emphasise the
  hovered side, say — throws the marker to the top-left corner of the map.
  Emphasise with colour, `box-shadow` or size instead.
- **Rebuilding a marker on its own hover eats the click.** `setHotSide` used to
  call `drawSideNums`, which removes and recreates every badge; the one under
  the pointer was destroyed between mousedown and mouseup, so the click never
  fired and the side never toggled. Repaint the class on the existing element
  via `getElement()` instead.
- **The middle of a side is already occupied.** The hollow "add a corner"
  handle sits exactly there, so a side badge dropped on the midpoint steals its
  clicks — or loses its own. The badge is nudged 20px outward from the shape's
  centre, in screen space, which is why `drawSideNums` runs again on `zoomend`.
- **`%TEMP%` does not exist off-Windows.** Use `tempfile.gettempdir()`.
- **`.sh` files need LF.** `.gitattributes` enforces it; without it they die on
  Linux with `bad interpreter`.

## Rebuilding this from scratch

The order that worked, if you ever need to reproduce it:

1. Verify the free imagery endpoints actually respond before designing around
   them — Esri World Imagery to z20 keyless, USGS via WMS (its `/tile/` REST
   path 404s), Nominatim for geocoding.
2. Write and validate section 1 against the closed form **before** any UI.
   Getting the maths right first is what makes everything after it trustworthy.
3. Leaflet + raster tiles, drawing tools hand-rolled. No Turf — the local-frame
   projection makes plain planar geometry correct at this scale.
4. Export/import next, tested by round-trip and by parsing a real Google Earth
   file (nested `Folder`, indented multi-line `<coordinates>`, `<tessellate>`).
5. Only then the desktop shell, and only because browser downloads are worse
   than native Save dialogs.

Verify each layer numerically before building the next one on top of it.
