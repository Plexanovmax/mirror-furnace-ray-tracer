#!/usr/bin/env python3
"""Interactive 2D and 3D ray tracing for PSI drawing DRW-00322318.
Run with Python 3, then open the printed local URL. No packages required.
To replace the overlay, keep furnace_drawing.png at its original crop.
"""

import json
import math
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from socketserver import TCPServer
from threading import Timer
from urllib.parse import parse_qs, urlsplit
import webbrowser


DEFAULT_RAY_COUNT = 12
RAY_COUNT_LIMITS = (2, 100)
MAX_REFLECTIONS = 12
LEVER_MM = 100.0
EPSILON = 1e-6
POSITION_LIMITS = {"x": (-2000, 2000), "y": (-2000, 2000), "z": (-2500, 2500)}
TILT_LIMITS = {"t": (-200, 200), "v": (-200, 200)}
LIMITS = {}
for group in ("s", "u", "l", "p"):
    dimensions = POSITION_LIMITS | TILT_LIMITS if group in ("u", "l") else POSITION_LIMITS
    LIMITS.update({group + axis: bounds for axis, bounds in dimensions.items()})
LIMITS.update(pz=(-4000, 4000), pr=(30, 550))
DEFAULTS = dict.fromkeys(LIMITS, 0)
DEFAULTS["pr"] = 50

OVERLAY_PATH = Path(__file__).resolve().with_name("furnace_drawing.png")
# PDF crop calibration at 1960 × 2650 px: X=0 at column 869.5,
# Z=0 at row 1511.5, with 1.2337 pixels per millimetre.
OVERLAY_BOUNDS = [-869.5 / 1.2337, 1511.5 / 1.2337, 1960 / 1.2337, 2650 / 1.2337]


@dataclass(frozen=True)
class Mirror:
    axial: float
    radial: float
    side: int
    rim_depth: float  # Distance inward from the ideal ellipse vertex.
    pole_depth: float

    @property
    def focal_distance(self):
        return math.sqrt(self.axial**2 - self.radial**2)

    @property
    def local_rim(self):
        return self.side * (self.axial - self.rim_depth)

    @property
    def nominal_rim(self):
        return self.side * self.focal_distance + self.local_rim


# Upper rim is 15 + 27.62 mm inward from the outer focus. Upper and
# lower physical axial heights are 179.7 and 184.8 mm respectively.
UPPER = Mirror(659.5, 416.4, 1, 148.08 + 15 + 27.62, 11.0)
LOWER = Mirror(470.5, 252.3, -1, 73.37 + 123.45, 73.37 + 123.45 - 184.8)
SOURCE_HEIGHTS = {"bottom": -2 * LOWER.focal_distance, "top": 2 * UPPER.focal_distance}
SOURCE_MODES = {"bottom": ("bottom",), "top": ("top",),
                "both": ("bottom", "top"), "none": ()}


class Pose:
    def __init__(self, mirror, x, y, z, tilt_xz, tilt_yz):
        self.mirror = mirror
        self.pivot = (x, y, mirror.nominal_rim + z)
        across = math.atan2(tilt_xz, LEVER_MM)
        depth = math.atan2(tilt_yz, LEVER_MM)
        c, s = math.cos(across), math.sin(across)
        d, t = math.cos(depth), math.sin(depth)
        # Apply the XZ rotation first, then the YZ rotation about the rim centre.
        self.rotation = ((c, 0, -s), (-t*s, d, -t*c), (d*s, t, d*c))
        self.inverse = tuple(zip(*self.rotation))

    @staticmethod
    def rotate(matrix, vector):
        return tuple(sum(a*b for a, b in zip(row, vector)) for row in matrix)

    def point(self, x, y, z):
        offset = self.rotate(self.rotation, (x, y, z - self.mirror.local_rim))
        return tuple(p + d for p, d in zip(self.pivot, offset))

    def profiles(self):
        m = self.mirror
        first = math.acos(1 - m.pole_depth / m.axial)
        last = math.acos(1 - m.rim_depth / m.axial)
        angles = [first + (last - first) * i / 79 for i in range(80)]
        profiles = [[self.point(sign * m.radial * math.sin(t), 0,
                                m.side * m.axial * math.cos(t))
                     for t in angles] for sign in (-1, 1)]
        if abs(self.rotation[1][2]) > EPSILON:
            for t in (first, last):
                radius = m.radial * math.sin(t)
                z = m.side * m.axial * math.cos(t)
                profiles.append([self.point(radius * math.cos(math.tau*i/96),
                                            radius * math.sin(math.tau*i/96), z)
                                 for i in range(97)])
        return profiles

    def inner_focus(self):
        return self.point(0, 0, -self.mirror.side * self.mirror.focal_distance)

    def intersection(self, origin, direction):
        m = self.mirror
        offset = tuple(p - pivot for p, pivot in zip(origin, self.pivot))
        px, py, pz = self.rotate(self.inverse, offset)
        pz += m.local_rim
        dx, dy, dz = self.rotate(self.inverse, direction)
        a = (dx*dx + dy*dy) / m.radial**2 + (dz / m.axial)**2
        b = (px*dx + py*dy) / m.radial**2 + pz * dz / m.axial**2
        c = (px*px + py*py) / m.radial**2 + (pz / m.axial)**2 - 1
        discriminant = b*b - a*c
        if discriminant < 0:
            return None

        root = math.sqrt(discriminant)
        for distance in ((-b - root) / a, (-b + root) / a):
            if distance <= EPSILON:
                continue
            hit_x = px + distance * dx
            hit_y = py + distance * dy
            hit_z = pz + distance * dz
            depth = m.axial - m.side * hit_z
            if not m.pole_depth - EPSILON <= depth <= m.rim_depth + EPSILON:
                continue

            # The ellipse gradient points outward. The concave coating
            # reflects only rays approaching from the interior.
            gradient = (hit_x / m.radial**2, hit_y / m.radial**2, hit_z / m.axial**2)
            length = math.hypot(*gradient)
            normal = self.rotate(self.rotation, tuple(n / length for n in gradient))
            dot = sum(v*n for v, n in zip(direction, normal))
            if dot <= EPSILON:
                return distance, None
            return distance, tuple(v - 2*dot*n for v, n in zip(direction, normal))
        return None


def emission_directions(spatial, ray_count=DEFAULT_RAY_COUNT):
    if not spatial:
        return [(math.cos(a), 0, math.sin(a))
                for a in (math.tau*i / ray_count for i in range(ray_count))]
    if ray_count == DEFAULT_RAY_COUNT:
        # Preserve the original twelve directions at the vertices of an icosahedron.
        phi = (1 + math.sqrt(5)) / 2
        length = math.hypot(1, phi)
        return [tuple(v / length for v in direction)
                for a in (-1, 1) for b in (-phi, phi)
                for direction in ((0, a, b), (a, b, 0), (b, 0, a))]
    if ray_count == 2:
        return [(0, 0, 1), (0, 0, -1)]

    # Equal-area latitude bands with golden-angle azimuth steps.
    golden_angle = math.pi * (3 - math.sqrt(5))
    directions = []
    for index in range(ray_count):
        z = 1 - (2*index + 1) / ray_count
        radial = math.sqrt(1 - z*z)
        angle = index * golden_angle
        directions.append((radial * math.cos(angle), radial * math.sin(angle), z))
    return directions


def ray_crossings(rays, height):
    crossings = []
    for number, ray in enumerate(rays, 1):
        points, tail = ray["points"], ray["tail"]
        relative = [point[2] - height for point in points]

        def add(x, y, vertical):
            crossings.append(dict(x=x, y=y, up=vertical > 0, ray=number))

        for i, (start_z, end_z) in enumerate(zip(relative, relative[1:])):
            if min(start_z, end_z) < -EPSILON and max(start_z, end_z) > EPSILON:
                start, end = points[i:i+2]
                fraction = -start_z / (end_z - start_z)
                add(*(a + fraction*(b - a) for a, b in zip(start[:2], end[:2])),
                    end_z - start_z)
        for i, z in enumerate(relative[1:], 1):
            if abs(z) > EPSILON:
                continue
            # A reflection on the plane counts only if the ray crosses it.
            after = relative[i+1] if i+1 < len(relative) else tail[2] if tail else 0
            before = relative[i-1]
            if ((before < -EPSILON and after > EPSILON)
                    or (before > EPSILON and after < -EPSILON)):
                add(*points[i][:2], after)
        if tail and abs(tail[2]) > 1e-12:
            distance = -relative[-1] / tail[2]
            if abs(relative[-1]) > EPSILON and distance > 0:
                add(points[-1][0] + distance*tail[0], points[-1][1] + distance*tail[1], tail[2])
    return crossings


def simulate(values, spatial=False, ray_count=DEFAULT_RAY_COUNT, *,
             source_mode="bottom", mirrors_enabled=(True, True)):
    upper, lower = [Pose(mirror, values[p + "x"], values[p + "y"] if spatial else 0,
                         values[p + "z"], values[p + "t"], values[p + "v"] if spatial else 0)
                    for mirror, p in ((UPPER, "u"), (LOWER, "l"))]
    sources = {location: (values["sx"], values["sy"] if spatial else 0,
                          SOURCE_HEIGHTS[location] + values["sz"])
               for location in SOURCE_MODES[source_mode]}
    poses = (upper, lower)
    active_mirrors = [(i, pose) for i, pose in enumerate(poses) if mirrors_enabled[i]]
    rays = []
    per_source, remainder = divmod(ray_count, max(1, len(sources)))
    emissions = ((source, direction) for index, source in enumerate(sources.values())
                 for direction in emission_directions(spatial, per_source + (index < remainder)))
    for source, direction in emissions:
        points = [source]
        reflections = []
        tail = None
        end = "limit"
        for _ in range(MAX_REFLECTIONS):
            hits = [(hit[0], number, hit[1]) for number, pose in active_mirrors
                    if (hit := pose.intersection(points[-1], direction)) is not None]
            if not hits:
                tail = direction
                end = "escaped"
                break
            distance, number, reflected = min(hits)
            points.append(tuple(p + distance*v for p, v in zip(points[-1], direction)))
            if reflected is None:
                end = "back"
                break
            reflections.append(number)
            direction = reflected
        rays.append(dict(points=points, tail=tail, reflections=reflections, end=end))
    # The furnace display is the XZ projection. Detector intersections are
    # evaluated on the full 3D paths before projecting them.
    detector = None
    if spatial:
        px, py, pz = (values[axis] for axis in ("px", "py", "pz"))
        radius = values.get("pr", DEFAULTS["pr"])
        plane_crossings = ray_crossings(rays, pz)
        crossings = [hit for hit in plane_crossings
                     if math.hypot(hit["x"] - px, hit["y"] - py) <= radius + EPSILON]
        detector = dict(centre=[px, py, pz], radius=radius, crossings=crossings,
                        missed=len(plane_crossings) - len(crossings))
    def project(point):
        return point[0], point[2]

    for ray in rays:
        ray["points"] = [project(point) for point in ray["points"]]
        if ray["tail"] is not None:
            ray["tail"] = project(ray["tail"])
    projected_sources = {location: project(source) for location, source in sources.items()}
    return dict(sources=projected_sources, source=next(iter(projected_sources.values()), None),
                mirrors=[[list(map(project, profile)) for profile in pose.profiles()] if enabled else []
                         for pose, enabled in zip(poses, mirrors_enabled)],
                foci=[project(pose.inner_focus()) if enabled else None
                      for pose, enabled in zip(poses, mirrors_enabled)],
                rays=rays, detector=detector)


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Mirror furnace ray tracer</title>
<style>
:root { color-scheme: light; font:15px/1.45 system-ui,sans-serif; color:#193142; background:#eff4f6; }
* { box-sizing:border-box; }
body { margin:0; }
header { padding:24px 30px 18px; border-bottom:1px solid #d5e0e6; background:#fff; }
h1 { font-size:clamp(23px,3vw,31px); letter-spacing:-.035em; margin:0; }
header p, section p { margin:3px 0 0; color:#496375; font-size:13px; }
main { display:grid; grid-template-columns:minmax(0,1fr) 365px; gap:16px; padding:16px;
       max-width:1850px; margin:auto; align-items:start; }
.views { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1.8fr) minmax(0,1fr);
    grid-template-areas:"zoom furnace furnace";
    min-width:0; gap:12px; align-items:start; }
.views:has(#detector:not([hidden])) { grid-template-areas:"zoom furnace detector"; }
#zoom-panel { grid-area:zoom; }
#furnace-panel { grid-area:furnace; }
#detector { grid-area:detector; }
[hidden] { display:none !important; }
.plot, aside { border:1px solid #d5e0e6; border-radius:14px; background:#fff; overflow:hidden; }
.toolbar, .summary { display:flex; align-items:center; justify-content:space-between;
    flex-wrap:wrap; gap:8px 12px; padding:11px 14px; }
.toolbar { border-bottom:1px solid #e1e9ed; }
.summary { border-top:1px solid #e1e9ed; color:#496375; font-size:13px; justify-content:flex-start; }
.summary span { white-space:nowrap; }
.key { display:inline-block; width:16px; border-top:3px solid var(--swatch); margin-right:5px; vertical-align:middle; }
#plot, #spots, #zoom { display:block; width:100%; background:#f8fbfc; }
#plot { height:clamp(490px,72vh,1000px); cursor:crosshair; }
#zoom, #spots { aspect-ratio:1; }
.furnace-canvas { position:relative; overflow:hidden; }
#lens { position:absolute; pointer-events:none; transform:translate(-50%,-50%);
    border:1.5px solid #176c75; box-shadow:0 0 0 1px #ffffffc9; background:#59c6b814; }
.range-row, .toolbar .range-row { display:grid; grid-template-columns:auto minmax(50px,1fr) auto;
    gap:8px; width:100%; align-items:center; margin-top:5px; }
.range-row output { min-width:3ch; text-align:right; font-variant-numeric:tabular-nums; }
.toolbar label { display:inline-flex; align-items:center; gap:6px; }
#opacity { width:88px; }
#opacity-value { min-width:35px; text-align:right; font-variant-numeric:tabular-nums; }
.check-label { display:flex; align-items:center; gap:7px; margin:0; }
.check-label input { margin:0; accent-color:#126f78; }
.controls label.inactive, .angle.inactive { color:#778892; }
input:disabled { opacity:.45; }
section { padding:12px 14px 14px; border-bottom:1px solid #e1e9ed; }
section:last-child { border:0; }
h2 { font-size:15px; margin:0; }
.controls { display:grid; grid-template-columns:61px minmax(55px,1fr) 30px 63px 30px;
    align-items:center; gap:7px 6px; margin-top:8px; }
label { font-size:13px; }
button, input, select { font:inherit; }
button, select, input[type=number] { border:1px solid #b9cbd5; color:#193142; background:#fff;
    border-radius:6px; min-height:31px; }
button, select { padding:3px 8px; cursor:pointer; }
button:hover { background:#e9f4f5; }
button:disabled { opacity:.38; cursor:default; }
input[type=range] { accent-color:#126f78; margin:0; width:100%; min-width:0; }
input[type=number] { width:100%; padding:3px 5px; text-align:right; appearance:textfield;
    font-variant-numeric:tabular-nums; }
input::-webkit-inner-spin-button { appearance:none; }
.angle { grid-column:2/-1; text-align:right; color:#496375; font-size:12px; }
.note { margin:0; padding:11px 14px; color:#496375; font-size:12px; }
@media (max-width:1200px) {
    main { grid-template-columns:minmax(0,1fr); padding:10px; }
    header { padding:15px; }
    #plot { height:clamp(405px,64vh,650px); }
}
@media (max-width:760px) {
    .views { grid-template-columns:minmax(0,1fr); grid-template-areas:"furnace" "zoom"; }
    .views:has(#detector:not([hidden])) { grid-template-areas:"furnace" "zoom" "detector"; }
}
</style>
</head>
<body>
<header>
  <h1>Mirror furnace ray tracer</h1>
  <p id="mode">DRW-00322318 · XZ cross-section</p>
</header>
<main>
  <div class="views">
  <div class="plot" id="furnace-panel">
    <div class="toolbar">
      <strong>Optical paths</strong>
      <label for="opacity">PDF opacity <input id="opacity" type="range" min="0" max="100" value="30" step="1">
        <output id="opacity-value" for="opacity">30%</output></label>
      <label for="view">View
      <select id="view"><option value="fit">Full furnace</option>
        <option value="focus">Focus detail</option>
        <option value="bottom">Bottom source detail</option>
        <option value="top">Top source detail</option>
        <option value="plane">Plane height detail</option></select></label>
    </div>
    <div class="furnace-canvas">
      <canvas id="plot" aria-label="Furnace side view with mirrors and ray paths"></canvas>
      <div id="lens" aria-hidden="true"></div>
    </div>
    <div class="summary">
      <span id="key-u"><i class="key" style="--swatch:#126b8f"></i>Upper mirror</span>
      <span id="key-l"><i class="key" style="--swatch:#704691"></i>Lower mirror</span>
      <span id="key-s"><i class="key" style="--swatch:#bd4917"></i>Rays</span>
      <span id="plane-key" hidden><i class="key" style="--swatch:#287870"></i>Detector edge</span>
      <span id="status" role="status">Tracing…</span>
    </div>
    <p class="note">The dashed extensions show escaping rays. The dotted ellipse outlines show the nominal geometry.
        The cross marks the nominal common focus. Rings mark the current mirror foci.
        Shutters and mechanical parts are omitted.</p>
  </div>
  <div class="plot" id="zoom-panel">
    <div class="toolbar">
      <strong>Hover zoom</strong>
      <label class="range-row" for="lens-size">Lens size
        <input id="lens-size" type="range" min="20" max="180" value="70" step="1">
        <output id="lens-size-value" for="lens-size">70 px</output>
      </label>
    </div>
    <canvas id="zoom" aria-label="Magnified view of the region under the furnace lens"></canvas>
    <div class="summary"><span id="zoom-scale"></span></div>
    <p class="note">Hover over the furnace to move the lens. A smaller lens gives more magnification.</p>
  </div>
  <div class="plot" id="detector" hidden>
    <div class="toolbar">
      <strong>Plane top view</strong>
      <label for="spot-radius">View
        <select id="spot-radius">
          <option value="fit">Fit detector</option>
          <option value="20">±20 mm</option><option value="50">±50 mm</option>
          <option value="100">±100 mm</option><option value="250">±250 mm</option>
          <option value="500">±500 mm</option><option value="1000">±1000 mm</option>
          <option value="2000">±2000 mm</option>
        </select>
      </label>
    </div>
    <canvas id="spots" aria-label="XY plane top view with ray crossings and 10 mm rings"></canvas>
    <div class="summary">
      <span><i class="key" style="--swatch:#bd4917"></i>Downward</span>
      <span><i class="key" style="--swatch:#126b8f"></i>Upward</span>
      <span id="spot-status" role="status"></span>
    </div>
    <p class="note">Rings are spaced every 10 mm. Counts include repeated crossings and coincident spots.</p>
  </div>
  </div>
  <aside>
    <section><label class="check-label"><input type="checkbox" id="fine-steps">Fine steps (0.1 mm)</label></section>
    <section><h2><label class="check-label"><input type="checkbox" id="plane-toggle">XY detector plane · 3D</label></h2>
      <p id="plane-help">Enable to trace in 3D and unlock the Y controls.</p><div class="controls" id="p" hidden></div></section>
    <section><h2>Point source</h2>
      <label class="range-row" for="source-mode">Location
        <select id="source-mode">
          <option value="bottom" selected>Bottom</option>
          <option value="top">Top</option>
          <option value="both">Both</option>
          <option value="none">None</option>
        </select>
      </label>
      <label class="range-row" for="ray-count">Rays
        <input id="ray-count" type="range" min="__RAY_MIN__" max="__RAY_MAX__" value="__RAY_DEFAULT__" step="1">
        <output id="ray-count-value" for="ray-count">__RAY_DEFAULT__</output>
      </label>
      <p>Offsets apply to each active source from its mirror’s nominal outer focus. Both splits the ray count between the sources.</p>
      <div class="controls" id="s"></div></section>
    <section><h2><label class="check-label"><input type="checkbox" id="enable-u" checked>Upper mirror</label></h2>
      <p>Offsets from its nominal rim centre.</p><div class="controls" id="u"></div></section>
    <section><h2><label class="check-label"><input type="checkbox" id="enable-l" checked>Lower mirror</label></h2>
      <p>Offsets from its nominal rim centre.</p><div class="controls" id="l"></div></section>
    <section><button id="reset" type="button">Reset to drawing</button>
      <p>Positive X is right and positive Z is up in the furnace view. Tilts rotate about the shifted rim centre.
          Enter rise on a 100 mm lever. XZ tilt is applied before YZ tilt.</p></section>
  </aside>
</main>
<script>
const limits = __LIMITS__;
const defaults = __DEFAULTS__;
const state = {...defaults};
const controls = {};
const nominal = __NOMINAL__;
const overlayBounds = __OVERLAY_BOUNDS__;
const names = {s:'Point source', u:'Upper mirror', l:'Lower mirror', p:'Plane'};
const objectToggles = ['u', 'l'].map(key => [key, document.querySelector('#enable-' + key)]);
const enabled = Object.fromEntries(objectToggles);
const colors = ['#126b8f', '#704691'];
const plot = document.querySelector('#plot');
const lens = document.querySelector('#lens');
const lensSize = document.querySelector('#lens-size');
const lensSizeValue = document.querySelector('#lens-size-value');
const zoom = document.querySelector('#zoom');
const zoomScale = document.querySelector('#zoom-scale');
const rayCount = document.querySelector('#ray-count');
const rayCountValue = document.querySelector('#ray-count-value');
const status = document.querySelector('#status');
const view = document.querySelector('#view');
const sourceMode = document.querySelector('#source-mode');
const sourceDetails = ['bottom', 'top'].map(location =>
    [location, view.querySelector(`option[value=${location}]`)]);
const fineSteps = document.querySelector('#fine-steps');
const planeToggle = document.querySelector('#plane-toggle');
const planeControls = document.querySelector('#p');
const planeHelp = document.querySelector('#plane-help');
const planeDetail = view.querySelector('option[value=plane]');
const detectorPanel = document.querySelector('#detector');
const spots = document.querySelector('#spots');
const spotRadius = document.querySelector('#spot-radius');
const opacity = document.querySelector('#opacity');
const opacityValue = document.querySelector('#opacity-value');
const overlay = new Image();
let geometry = nominal;
let requestVersion = 0;
let scheduled = false;
let lensPointer = null;
let mainFrame;
let zoomScheduled = false;

const movementStep = () => fineSteps.checked ? 0.1 : 1;

function refreshControl(key) {
    const control = controls[key];
    const {slider, field, minus, plus, angle, label, name} = control;
    const value = state[key];
    // Keep a convenient local slider span while allowing larger typed offsets.
    control.span = Math.max(control.span, Math.ceil(Math.abs(value) / 50) * 50);
    slider.min = Math.max(limits[key][0], -control.span);
    slider.max = Math.min(limits[key][1], control.span);
    slider.value = field.value = value;
    const disabled = (key.startsWith('s') && sourceMode.value === 'none')
        || enabled[key[0]]?.checked === false
        || (!planeToggle.checked && (key.endsWith('y') || key.endsWith('v') || key.startsWith('p')));
    slider.disabled = field.disabled = disabled;
    label.classList.toggle('inactive', disabled);
    minus.disabled = disabled || value <= limits[key][0];
    plus.disabled = disabled || value >= limits[key][1];
    minus.setAttribute('aria-label', `Decrease ${name} by ${movementStep()} mm`);
    plus.setAttribute('aria-label', `Increase ${name} by ${movementStep()} mm`);
    if (angle) {
        angle.classList.toggle('inactive', disabled);
        const degrees = Math.atan2(value, 100) * 180 / Math.PI;
        angle.textContent = `Rotation ${degrees.toFixed(2)}°`;
    }
}

function update(key, value) {
    const [minimum, maximum] = limits[key];
    const next = Math.max(minimum, Math.min(maximum, Math.round(value * 10) / 10));
    const changed = state[key] !== next;
    state[key] = next;
    refreshControl(key);
    if (changed) queueTrace();
}

for (const [key, [minimum, maximum]] of Object.entries(limits)) {
    const group = document.querySelector('#' + key[0]);
    const component = {x:'X', y:'Y', z:'Z', t:'Tilt XZ', v:'Tilt YZ', r:'Radius'}[key[1]];
    const name = names[key[0]] + ' ' + component;
    const label = document.createElement('label');
    label.htmlFor = key;
    label.textContent = component + ' [mm]';

    const slider = Object.assign(document.createElement('input'), {type:'range', step:1});
    slider.setAttribute('aria-label', name + ' slider (mm)');
    slider.addEventListener('input', () => update(key, slider.valueAsNumber));

    const field = Object.assign(document.createElement('input'),
        {type:'number', id:key, min:minimum, max:maximum, step:0.1});
    field.setAttribute('aria-label', name + ' (mm)');
    function commit() {
        if (Number.isFinite(field.valueAsNumber)) update(key, field.valueAsNumber);
        else refreshControl(key);
    }
    field.addEventListener('change', commit);
    field.addEventListener('keydown', event => {
        if (event.key === 'Enter') commit();
        if (event.key === 'Escape') refreshControl(key);
    });

    const plus = Object.assign(document.createElement('button'), {type:'button', textContent:'+'});
    const minus = Object.assign(document.createElement('button'), {type:'button', textContent:'−'});
    plus.addEventListener('click', () => update(key, state[key] + movementStep()));
    minus.addEventListener('click', () => update(key, state[key] - movementStep()));
    group.append(label, slider, minus, field, plus);
    const angle = 'tv'.includes(key[1]) ? document.createElement('div') : null;
    if (angle) {
        angle.className = 'angle';
        group.append(angle);
    }
    const initialSpan = key === 'pz' || key === 'pr' ? maximum : 150;
    controls[key] = {name, label, slider, field, plus, minus, angle, initialSpan, span:initialSpan};
    refreshControl(key);
}

function queueTrace() {
    requestVersion += 1;
    if (scheduled) return;
    scheduled = true;
    requestAnimationFrame(async () => {
        scheduled = false;
        const version = requestVersion;
        try {
            const query = new URLSearchParams({
                ...state, rays:rayCount.value, spatial:Number(planeToggle.checked), source_mode:sourceMode.value
            });
            for (const [key, checkbox] of objectToggles) query.set('enable_' + key, Number(checkbox.checked));
            const response = await fetch('/trace?' + query, {cache:'no-store'});
            if (!response.ok) throw new Error(response.status);
            const result = await response.json();
            if (version === requestVersion) {
                geometry = result;
                draw();
            }
        } catch {
            if (version === requestVersion) status.textContent = 'Connection lost. Check the Python server.';
        }
    });
}

function frame(width, height) {
    let centreX, centreZ, scale;
    const detector = planeToggle.checked ? geometry.detector : null;
    if (view.value === 'fit') {
        const points = [...Object.values(geometry.sources), ...geometry.foci, ...geometry.mirrors.flat(2),
            ...nominal.guides.flat(), ...geometry.rays.flatMap(ray => ray.points)].filter(Boolean);
        if (detector) {
            const [px, , pz] = detector.centre;
            points.push([px - detector.radius,pz], [px + detector.radius,pz]);
        }
        const xs = points.map(p => p[0]), zs = points.map(p => p[1]);
        const minX = Math.min(...xs), maxX = Math.max(...xs);
        const minZ = Math.min(...zs), maxZ = Math.max(...zs);
        centreX = (minX + maxX) / 2;
        centreZ = (minZ + maxZ) / 2;
        scale = Math.min(width / (maxX - minX + 190), height / (maxZ - minZ + 240));
    } else {
        if (geometry.sources[view.value]) [centreX, centreZ] = geometry.sources[view.value];
        else if (view.value === 'plane' && detector) [centreX, , centreZ] = detector.centre;
        else [centreX, centreZ] = [0, 0];
        const span = view.value === 'plane' && detector ? Math.max(60, 2.4*detector.radius) : 60;
        scale = Math.min(width / span, height / 60);
    }
    return makeFrame(width, height, centreX, centreZ, scale);
}

function makeFrame(width, height, centreX, centreZ, scale) {
    return {scale, x: x => width / 2 + (x - centreX) * scale,
        z: z => height / 2 - (z - centreZ) * scale,
        bounds:[centreX-width/(2*scale), centreX+width/(2*scale),
                centreZ-height/(2*scale), centreZ+height/(2*scale)]};
}

function prepareCanvas(canvas) {
    const width = canvas.clientWidth, height = canvas.clientHeight;
    const density = window.devicePixelRatio || 1;
    canvas.width = Math.round(width * density);
    canvas.height = Math.round(height * density);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(density, 0, 0, density, 0, 0);
    return {width, height, ctx};
}

function drawFurnace({width, height, ctx:context}, screen) {
    const {x, z, scale, bounds:[xmin, xmax, zmin, zmax]} = screen;
    if (overlay.complete && overlay.naturalWidth && opacity.valueAsNumber > 0) {
        const [left, top, imageWidth, imageHeight] = overlayBounds;
        context.save();
        context.globalAlpha = opacity.valueAsNumber / 100;
        context.drawImage(overlay, x(left), z(top), scale * imageWidth, scale * imageHeight);
        context.restore();
    }
    function path(points, color, lineWidth, dash = []) {
        context.beginPath();
        points.forEach(([px,pz], index) => index ? context.lineTo(x(px), z(pz)) : context.moveTo(x(px), z(pz)));
        context.lineJoin = 'round';
        context.strokeStyle = color;
        context.lineWidth = lineWidth;
        context.setLineDash(dash);
        context.stroke();
        context.setLineDash([]);
    }
    const steps = [0.25, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 250, 500, 1000, 2000];
    const spacing = steps.find(step => step * scale >= 47) || steps.at(-1);
    context.font = '11px system-ui';
    context.fillStyle = '#4d6575';
    context.textAlign = 'center';
    for (let value = Math.ceil(xmin/spacing)*spacing; value <= xmax; value += spacing) {
        path([[value,zmin],[value,zmax]], value === 0 ? '#b5c8d1' : '#e1e9ed', 1);
        if (x(value) > 44 && x(value) < width - 59) context.fillText(value, x(value), height - 10);
    }
    context.textAlign = 'left';
    for (let value = Math.ceil(zmin/spacing)*spacing; value <= zmax; value += spacing) {
        path([[xmin,value],[xmax,value]], value === 0 ? '#b5c8d1' : '#e1e9ed', 1);
        if (z(value) > 25 && z(value) < height - 30) context.fillText(value, 8, z(value) - 4);
    }
    context.fillText('Z [mm]', 8, 16);
    context.fillText('X [mm]', width - 42, height - 10);

    const detector = planeToggle.checked ? geometry.detector : null;
    if (detector) {
        const [px, , pz] = detector.centre;
        const left = px - detector.radius, right = px + detector.radius;
        path([[left,pz],[right,pz]], '#59c6b872', 10);
        path([[left,pz],[right,pz]], '#216c64', 2);
        for (const edge of [left, right]) {
            path([[edge,pz-5/scale],[edge,pz+5/scale]], '#216c64', 1.6);
        }
    }
    nominal.guides.forEach((guide, index) => {
        if (geometry.mirrors[index].length) path(guide, '#aabcc5', 1, [2,5]);
    });
    for (const ray of geometry.rays) {
        path(ray.points, '#bd4917b8', 1.45);
        const projectedLength = ray.tail ? Math.hypot(...ray.tail) : 0;
        if (projectedLength > 1e-12) {
            const [tx,tz] = ray.points.at(-1);
            const tailLength = (1 + Math.hypot(Math.max(Math.abs(tx-xmin), Math.abs(tx-xmax)),
                Math.max(Math.abs(tz-zmin), Math.abs(tz-zmax)))) / projectedLength;
            path([[tx,tz],[tx+tailLength*ray.tail[0],tz+tailLength*ray.tail[1]]],
                 '#bd491789', 1.2, [5,5]);
        }
    }
    geometry.mirrors.forEach((mirror, index) => mirror.forEach(profile => path(profile, colors[index], 3.5)));

    function marker(px, pz, radius, color, filled) {
        context.beginPath();
        context.arc(x(px), z(pz), radius, 0, Math.PI * 2);
        context.strokeStyle = color;
        context.fillStyle = filled ? color : '#f8fbfc';
        context.lineWidth = 1.5;
        context.fill();
        context.stroke();
    }
    path([[-7/scale,0],[7/scale,0]], '#193142', 1.3);
    path([[0,-7/scale],[0,7/scale]], '#193142', 1.3);
    geometry.foci.forEach((focus, index) => {
        if (focus) marker(...focus, 4 + index * 3, colors[index], false);
    });
    Object.values(geometry.sources).forEach(source => marker(...source, 4.5, '#a53c12', true));
    if (detector) {
        const [px, , pz] = detector.centre;
        marker(px, pz, 3, '#216c64', true);
        if (z(pz) >= 0 && z(pz) <= height
                && x(px + detector.radius) >= 0 && x(px - detector.radius) <= width) {
            const label = `Detector · Z ${pz} mm · Ø ${2*detector.radius} mm`;
            const textWidth = context.measureText(label).width;
            const left = Math.max(5, Math.min(width - textWidth - 13, x(px) - textWidth/2));
            const top = z(pz) > 24 ? z(pz) - 25 : z(pz) + 9;
            context.fillStyle = '#f8fbfce8';
            context.fillRect(left, top, textWidth + 8, 16);
            context.fillStyle = '#16584f';
            context.textAlign = 'left';
            context.fillText(label, left + 4, top + 12);
        }
    }
}

function drawZoom() {
    const canvas = prepareCanvas(zoom);
    const width = plot.clientWidth, height = plot.clientHeight;
    if (!mainFrame || !width || !height || !canvas.width || !canvas.height) return;
    const {x, z, bounds, scale} = mainFrame;
    // Start at the common focus when visible. Keep the last hovered point
    // when the pointer leaves so the lens slider can adjust that region.
    const focusVisible = x(0) >= 0 && x(0) <= width && z(0) >= 0 && z(0) <= height;
    const left = lensPointer ? lensPointer[0]*width : focusVisible ? x(0) : width/2;
    const top = lensPointer ? lensPointer[1]*height : focusVisible ? z(0) : height/2;
    const pixels = lensSize.valueAsNumber;
    const side = pixels / scale;
    const zoomFactor = Math.min(canvas.width, canvas.height) / pixels;
    const centreX = bounds[0] + left/scale;
    const centreZ = bounds[3] - top/scale;
    lens.style.left = `${left}px`;
    lens.style.top = `${top}px`;
    lens.style.width = lens.style.height = `${pixels}px`;
    lensSizeValue.textContent = `${pixels} px`;
    zoomScale.textContent = `${zoomFactor.toFixed(1)}× · ${Number(side.toPrecision(3))} mm across lens`;
    drawFurnace(canvas, makeFrame(canvas.width, canvas.height, centreX, centreZ, scale*zoomFactor));
}

function moveLens(event) {
    if (event.pointerType === 'touch' && event.type === 'pointermove') return;
    const box = plot.getBoundingClientRect();
    lensPointer = [(event.clientX - box.left)/box.width, (event.clientY - box.top)/box.height];
    if (zoomScheduled) return;
    zoomScheduled = true;
    requestAnimationFrame(() => {
        zoomScheduled = false;
        drawZoom();
    });
}

function draw() {
    const canvas = prepareCanvas(plot);
    mainFrame = frame(canvas.width, canvas.height);
    drawFurnace(canvas, mainFrame);
    const count = geometry.rays.length;
    const reflected = geometry.rays.filter(ray => ray.reflections.length).length;
    const reachedBoth = geometry.rays.filter(ray => new Set(ray.reflections).size === 2).length;
    status.textContent = count
        ? `${reflected}/${count} reflected · ${reachedBoth}/${count} hit both mirrors`
        : 'Source off';
    drawZoom();
    drawDetector();
}

function drawDetector() {
    if (detectorPanel.hidden) return;
    const status = document.querySelector('#spot-status');
    const {width, height, ctx} = prepareCanvas(spots);
    if (!geometry.detector) {
        status.textContent = 'Tracing…';
        return;
    }
    const {centre:[px, py], radius:detectorRadius, crossings, missed} = geometry.detector;
    const radius = spotRadius.value === 'fit' ? detectorRadius : Number(spotRadius.value);
    const cx = width/2, cy = height/2;
    const scale = Math.max(1, Math.min(cx,cy) - 32) / radius;
    const extent = radius * scale;
    const rings = Math.floor(Math.min(radius, detectorRadius) / 10);
    const labelEvery = Math.max(1, Math.ceil(34 / (10 * scale)));

    ctx.font = '11px system-ui';
    ctx.fillStyle = '#496375';
    ctx.textAlign = 'right';
    if (detectorRadius <= radius) {
        ctx.beginPath();
        ctx.arc(cx, cy, detectorRadius * scale, 0, 2 * Math.PI);
        ctx.fillStyle = '#59c6b814';
        ctx.fill();
        ctx.strokeStyle = '#216c64';
        ctx.lineWidth = 1.6;
        ctx.stroke();
        ctx.fillStyle = '#496375';
    }
    for (let ring = 1; ring <= rings; ring++) {
        const major = ring === rings || ring % labelEvery === 0;
        if (ring * 10 !== detectorRadius) {
            ctx.beginPath();
            ctx.arc(cx, cy, ring * 10 * scale, 0, 2 * Math.PI);
            ctx.strokeStyle = major ? '#aabfc9' : '#dce7ec';
            ctx.lineWidth = major ? 1 : 0.7;
            ctx.stroke();
        }
        if (ring === rings || (major && rings - ring >= labelEvery)) {
            ctx.fillText(`${ring*10} mm`, cx - 5, cy - ring*10*scale - 3);
        }
    }
    ctx.strokeStyle = '#9db5bf';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(cx-extent, cy);
    ctx.lineTo(cx+extent, cy);
    ctx.moveTo(cx, cy-extent);
    ctx.lineTo(cx, cy+extent);
    ctx.stroke();
    ctx.fillText('+X', cx + extent + 21, cy + 4);
    ctx.fillText('+Y', cx + 15, cy - extent - 14);

    const visible = crossings.filter(hit => Math.hypot(hit.x - px, hit.y - py) <= radius + 1e-6);
    for (const up of [false, true]) {
        for (const hit of visible.filter(hit => hit.up === up)) {
            ctx.beginPath();
            ctx.arc(cx + (hit.x - px)*scale, cy - (hit.y - py)*scale, up ? 5 : 3, 0, 2*Math.PI);
            if (up) {
                ctx.strokeStyle = '#126b8f';
                ctx.lineWidth = 1.5;
                ctx.stroke();
            } else {
                ctx.fillStyle = '#bd4917';
                ctx.fill();
            }
        }
    }
    const outside = crossings.length - visible.length;
    const rays = new Set(crossings.map(hit => hit.ray)).size;
    status.textContent = geometry.rays.length
        ? `${crossings.length} hits · ${rays}/${geometry.rays.length} rays`
            + (missed ? ` · ${missed} missed disk` : '')
            + (outside ? ` · ${outside} outside view` : '')
        : '0 hits · source off';
}

function updateMode() {
    const active = planeToggle.checked;
    planeControls.hidden = detectorPanel.hidden = !active;
    document.querySelector('#plane-key').hidden = !active;
    planeDetail.disabled = !active;
    if (!active && view.value === 'plane') view.value = 'fit';
    planeHelp.textContent = active
        ? 'X and Y set the disk centre. Z is measured from the nominal common focus.'
        : 'Enable to trace in 3D and unlock the Y controls.';
    document.querySelector('#mode').textContent = active
        ? 'DRW-00322318 · XZ projection · Rays sampled across the sphere'
        : 'DRW-00322318 · XZ cross-section · Rays evenly spaced in angle';
    Object.keys(controls).forEach(refreshControl);
}

function updateObjects() {
    for (const [key, checkbox] of objectToggles) {
        document.querySelector('#key-' + key).hidden = !checkbox.checked;
    }
    for (const [location, option] of sourceDetails) {
        option.disabled = sourceMode.value !== location && sourceMode.value !== 'both';
        if (option.disabled && view.value === location) view.value = 'fit';
    }
    rayCount.disabled = sourceMode.value === 'none';
    document.querySelector('#key-s').hidden = rayCount.disabled;
    Object.keys(controls).forEach(refreshControl);
}

document.querySelector('#reset').addEventListener('click', () => {
    sourceMode.value = 'bottom';
    for (const [, checkbox] of objectToggles) checkbox.checked = checkbox.defaultChecked;
    for (const key of Object.keys(state)) {
        state[key] = defaults[key];
        controls[key].span = controls[key].initialSpan;
        refreshControl(key);
    }
    rayCount.value = rayCount.defaultValue;
    rayCountValue.textContent = rayCount.value;
    updateObjects();
    queueTrace();
});
for (const control of [sourceMode, ...Object.values(enabled)]) {
    control.addEventListener('change', () => {
        updateObjects();
        queueTrace();
    });
}
planeToggle.addEventListener('change', () => {
    geometry = nominal;
    updateMode();
    draw();
    queueTrace();
});
fineSteps.addEventListener('change', () => Object.keys(controls).forEach(refreshControl));
opacity.addEventListener('input', () => {
    opacityValue.textContent = `${opacity.value}%`;
    draw();
});
rayCount.addEventListener('input', () => {
    rayCountValue.textContent = rayCount.value;
    queueTrace();
});
lensSize.addEventListener('input', drawZoom);
plot.addEventListener('pointermove', moveLens);
plot.addEventListener('pointerdown', moveLens);
overlay.onload = draw;
overlay.onerror = () => {
    opacity.disabled = true;
    opacityValue.textContent = 'PNG missing';
};
overlay.src = '/furnace_drawing.png';
view.addEventListener('change', draw);
spotRadius.addEventListener('change', drawDetector);
const resize = new ResizeObserver(draw);
resize.observe(plot);
resize.observe(zoom);
resize.observe(spots);
updateMode();
updateObjects();
draw();
</script>
</body>
</html>
"""


def nominal_data():
    geometry = simulate(DEFAULTS)
    geometry["guides"] = [
        [(m.radial * math.sin(math.tau * i / 320),
          m.side * m.focal_distance + m.axial * math.cos(math.tau * i / 320))
         for i in range(321)]
        for m in (UPPER, LOWER)
    ]
    return geometry


class LocalServer(ThreadingHTTPServer):
    def server_bind(self):
        # The simulator only serves loopback traffic and needs no reverse DNS.
        TCPServer.server_bind(self)
        self.server_name, self.server_port = self.server_address


class Handler(BaseHTTPRequestHandler):
    def send_content(self, content, content_type):
        payload = content.encode("utf-8") if isinstance(content, str) else content
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self):
        target = urlsplit(self.path)
        if target.path == "/":
            self.send_content(self.server.page, "text/html; charset=utf-8")
        elif target.path == "/furnace_drawing.png":
            try:
                overlay = OVERLAY_PATH.read_bytes()
            except FileNotFoundError:
                self.send_error(404, "Missing furnace_drawing.png")
                return
            self.send_content(overlay, "image/png")
        elif target.path == "/trace":
            parameters = parse_qs(target.query)
            try:
                values = {key: float(parameters.get(key, [default])[0])
                          for key, default in DEFAULTS.items()}
                ray_count = int(parameters.get("rays", [DEFAULT_RAY_COUNT])[0])
            except ValueError:
                self.send_error(400, "Enter numeric coordinates and an integer ray count.")
                return
            if any(not LIMITS[key][0] <= value <= LIMITS[key][1]
                   for key, value in values.items()) or not RAY_COUNT_LIMITS[0] <= ray_count <= RAY_COUNT_LIMITS[1]:
                self.send_error(400, "Value outside control range.")
                return
            spatial = parameters.get("spatial", ["0"])[0] == "1"
            legacy_source = "bottom" if parameters.get("enable_s", ["1"])[0] == "1" else "none"
            source_mode = parameters.get("source_mode", [legacy_source])[0]
            if source_mode not in SOURCE_MODES:
                self.send_error(400, "Choose bottom, top, both or none.")
                return
            mirrors_enabled = tuple(
                parameters.get("enable_" + key, ["1"])[0] == "1" for key in ("u", "l")
            )
            result = simulate(values, spatial, ray_count, source_mode=source_mode,
                              mirrors_enabled=mirrors_enabled)
            self.send_content(json.dumps(result), "application/json")
        else:
            self.send_error(404)

    def log_message(self, *_):
        pass


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-browser", action="store_true", help="Print the URL without launching a browser.")
    args = parser.parse_args()
    with LocalServer(("127.0.0.1", 0), Handler) as server:
        server.page = PAGE.replace("__LIMITS__", json.dumps(LIMITS)).replace(
            "__DEFAULTS__", json.dumps(DEFAULTS)).replace(
            "__NOMINAL__", json.dumps(nominal_data())).replace(
            "__OVERLAY_BOUNDS__", json.dumps(OVERLAY_BOUNDS)).replace(
            "__RAY_MIN__", str(RAY_COUNT_LIMITS[0])).replace(
            "__RAY_MAX__", str(RAY_COUNT_LIMITS[1])).replace(
            "__RAY_DEFAULT__", str(DEFAULT_RAY_COUNT))
        url = f"http://127.0.0.1:{server.server_port}"
        print(f"Mirror furnace simulator: {url}\nPress Ctrl+C to stop.", flush=True)
        if not args.no_browser:
            Timer(0.25, webbrowser.open, args=(url,)).start()
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
