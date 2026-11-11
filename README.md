# Text Cube

**Live app: https://keithworrell.github.io/textcube/**

Type three letters and get one solid shape whose front, side and top views each show a different letter, like the *Gödel, Escher, Bach* cover. Export a watertight STL that goes straight into Bambu Studio or any other slicer.

## How it works

1. Each letter is drawn to a canvas and turned into a 2D signed distance field.
2. Each letter is extruded through the cube along its own axis (front → Z, side → X, top → Y). The solid is where all three extrusions overlap, combined with a smooth-max so the ridges where letters cross are filleted.
3. The surface is meshed with surface nets, then smoothed with Taubin smoothing, which polishes without shrinking the shape.
4. A shadow check projects the solid back along each axis and shows how much of each letter survived. Red marks the parts that were lost.

All of this runs in a Web Worker in the browser. There is no build server and no dependencies beyond three.js, which loads from a CDN.

## Controls

| Control | Trade-off |
|---|---|
| Font | Heavy fonts leave more material for the other two letters to cut into. |
| Fill the face | Stretches letters to the full square: more complete letters, less natural proportions. |
| Boldness | Sturdier and fills gaps, but closes counters (the holes in A, B, O…). |
| Letter corner rounding | Softens each letter's outline; too much looks melted. |
| Edge fillet | Rounds the ridges where letters cross; trims a little off the silhouettes. |
| Polish | Removes stair-steps; high values soften detail. |
| Resolution | Use Draft while exploring and Fine/Max before export. |
| One solid piece | Drops floating fragments so the print is a single body. |
| ↻ per letter | Rotates a letter 90°. Some combinations work much better rotated. |

## Printing

- The STL is in millimetres, Z-up, with the top letter facing up as printed.
- The model comes out about 1 mm under the nominal size because of the edge rounding.
- Expect to need supports under overhanging arms.

## Development

```
src/core.js        geometry engine (worker): SDFs, intersection, surface nets, smoothing, shadow check
src/template.html  UI, three.js preview, STL export
src/build.py       inlines core.js into the template → index.html
```

```
python3 src/build.py prod   # regenerates index.html (three.js from jsDelivr)
```

`index.html` is the built output that GitHub Pages serves. After editing anything in `src/`, rebuild and commit it.
