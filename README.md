# Warehouse Area Measure

Measure building footprints, roof areas and distances from satellite imagery.
Draw on the map, read the square footage, save the result as a `.kml` you can
open in Google Earth.

No account, no API key, no server. Nothing you draw leaves your machine.

## Get it running

**Easiest — download the app.** Grab `WarehouseAreaMeasure.exe` from the
[Releases](../../releases) page and double-click it. One file, its own window,
no Python needed. Windows only.

**From source, any OS.** No install step, no dependencies:

```bash
git clone <this repo>
cd warehouse-measure
python -m http.server 8123
```

Then open <http://localhost:8123>. On Windows you can double-click `start.bat`
instead, on Linux and macOS `./start.sh` — both do the same thing and open your
browser for you.

**Build the desktop app yourself** (gives you a native window and real Save/Open
dialogs instead of browser downloads):

```bash
pip install -r requirements.txt
python build_exe.py
```

That writes `dist/WarehouseAreaMeasure.exe`. It is unsigned, so Windows
SmartScreen will warn anyone who downloads it — "More info → Run anyway".

## Using it

1. Search the site address (or paste `32.7555, -97.3308`), and zoom in close.
2. **Area** — click each corner of the roof, then double-click or press `Enter`
   to close the shape. Hold `Shift` while clicking to snap each new wall square
   to the previous one; most warehouses are rectangular and this is faster and
   more accurate than eyeballing corners.
3. **Line** — click along a route for a distance. Every segment is measured as
   well as the total.
4. **Offset** — select a shape, set a distance (100 ft by default) and press
   Create. You get a new shape that distance out all the way around, with the
   ring area reported separately. Good for fire lanes, setbacks and laydown
   yards. It works on lines too, which gives you a corridor of a set width.
5. The `+` button next to an area flips it to **subtract**, so you can cut a
   courtyard or an excluded section out of the site total. Lines never count
   toward the area total.

**Editing:** drag a corner handle to move it, click a hollow midpoint handle to
add a corner, right-click a corner to delete it. Drag the blue grip at the
centre to slide the whole shape around — that is also the reliable way to grab
a line, which is too thin to hit accurately.

Everything autosaves to your browser, so closing the tab does not lose work.

### Shortcuts

| Key | |
|---|---|
| `A` / `L` | Area / Line tool |
| `Shift` | Hold while drawing to snap to right angles |
| `Enter` | Finish the current shape |
| `Esc` | Cancel drawing, or deselect |
| `Del` | Delete the selected shape |

## Files

**Save .kml** writes a Google Earth file: one Placemark per shape, named,
styled, with the measurements in the description. Areas become Polygons and
lines become LineStrings. Blue is included, red dashed is excluded, amber is a
line. **.geojson** is the same data for QGIS or anything else.

**Open file…** reads `.kml` and `.geojson` back, including files written by
Google Earth and QGIS — nested folders, multi-line coordinates, MultiPolygons
and MultiLineStrings all work. You can also drag a file onto the window.
Polygon holes (inner rings) are ignored; use a subtract shape instead.

**Print** produces a one-page sheet with the map and a table of every shape.
Print to PDF to send it on.

## Satellite imagery

Out of the box the map uses **Esri World Imagery** — free, no key, no account,
good to zoom 20 (about 6 in per pixel). Also included are **USGS imagery**
(US only, public domain) and **OpenStreetMap** for street context.

If you have a Google Maps Platform API key you can paste it under *Satellite
imagery* in the sidebar. It switches the map to Google's **Map Tiles API**,
which is sharper in some areas. Two things to know:

- Google bills per tile. Esri does not.
- If Google fails for any reason — bad key, quota, network — the map falls back
  to Esri on its own rather than going blank.

The key is stored in your browser's local storage and is sent to Google in tile
URLs, so restrict it to your own referrer in the Google Cloud console.

## Accuracy

Areas are computed in a local East/North metre frame built from the true WGS84
radii of curvature at each shape's own latitude, then measured with the
shoelace formula. Checked against the closed-form ellipsoidal area: **error is
below 0.0001%** at building scale.

For comparison, the spherical formula behind Google Maps' `computeArea` and
Turf's `area` — what most web measuring tools use — runs 0.67% high at the
equator, 0.12% high at 40°N and 0.41% low at 64°N. On a 400,000 sq ft warehouse
in Texas that is about 470 sq ft of error. This tool does not have it.

The offset tool places every vertex at exactly the requested distance; its area
matches `A + Pd + πd²` to within 0.02%, the residual being the polygon
approximation of the rounded corners.

The math is not the limiting factor. These are:

- **Building lean.** Satellite imagery is rarely straight overhead, so a tall
  warehouse leans and its roof outline sits offset from its footprint. Pick one
  and trace it consistently all the way around.
- **Imagery date.** A recent expansion may not appear. Cross-check the county
  parcel viewer when the number matters.
- **Your clicking.** At zoom 20 a pixel is about 6 inches.

The optional **roof pitch** field converts flat footprint to sloped roof
surface (`sqrt(1 + (pitch/12)^2)`). Warehouses are typically 0.25:12 to 0.5:12,
which adds well under 1%. The **unit rate** field multiplies the total by a
price per unit area for a quick estimate.

## Layout

```
index.html      markup and layout
app.css         styling, including the print stylesheet
app.js          geodesy, drawing tools, offset, KML/GeoJSON import and export
vendor/         Leaflet 1.9.4 (BSD-2-Clause, see vendor/LEAFLET-LICENSE.txt)
desktop.py      native-window shell: local server + OS file dialogs
build_exe.py    icon generation + PyInstaller build
start.bat       run it in a browser (Windows)
start.sh        run it in a browser (Linux / macOS)
```

`app.js` is organised in numbered sections; the geodesy is section 1 and is
independent of everything else.

## Licence

MIT — see [LICENSE](LICENSE). Leaflet is bundled under its own BSD-2-Clause
licence. Map imagery belongs to its providers and their attribution is
displayed on the map, as their terms require.
