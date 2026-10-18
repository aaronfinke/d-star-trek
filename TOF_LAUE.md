# Neutron TOF-Laue data

This version of d\*TREK processes neutron time-of-flight (TOF) Laue data of
the LANSCE PCS (cylindrical detector, the rest of this section) and of panel
instruments such as SNS MANDI (event data or NXlauetof histograms, see
[Panel instruments](#panel-instruments-sns-mandi-nxlauetof)).

## LANSCE PCS

This version of d\*TREK processes the wavelength-resolved time-of-flight (TOF)
Laue data of the Protein Crystallography Station (PCS) at LANSCE, as the lost
"heavily modified d\*TREK" of the PCS did (Langan & Greene, *J. Appl. Cryst.*
2004, **37**, 253–257). It was rebuilt from that paper and the outputs of the
original program, and it reads the original PCS files and the original
`.head` files. The command-line options of the original recipes
(`laue_process.txt`) are accepted. Everything is in the command-line
programs; dtdisplay and dtprocess are unchanged.

The TOF-Laue code paths run only for headers that describe TOF-Laue data
(`Detector1_ROTATION_ARRAY` and `Detector1_ROTATION_L1`). X-ray processing is
unchanged: dtfind and dtindex give identical results on X-ray images. The one
change that affects X-ray refinement is a bug fix: `dtrefine -CrysConstraints`
(an unconstrained cell for a non-triclinic space group) used the cell
derivatives of the constrained cell, so it converged slowly or to a wrong
cell, and printed the b and c values as their sigmas.

### Data files

A PCS run is one NeXus file. It holds eight detector segments of
96 TOF channels × 256 × 256 pixels, plus the TOF bin edges. The original
files are HDF4 (`*.nx.hdf`). d\*TREK reads HDF5, so convert them first:

    scripts/pcs_h4toh5.sh [-o outdir] R000466.nx.hdf ...

A `SCAN_TEMPLATE` or image name of `R000466.nx.hdf` then finds the
`R000466.h5` beside it. The reader stitches the segments into one
256 × 2048 cylindrical image per channel, the same layout as the original
program. TOF channel *c* is image *c* of the scan
(`SCAN_ROTATION=0 96 1`, `-seq 0 95`). The rotation fields of the reflection
lists hold TOF channels.

`dtpcsproject run.h5 out.img [-chan N | -sum [c0 c1] | -max [c0 c1]]` writes
one channel, the sum, or the maximum ("overlay") of a range of channels as an
ordinary image, for viewing.

### Processing one run

`scripts/pcs_process_run.sh` runs the steps of `laue_process.txt`:

    scripts/pcs_process_run.sh -o R000466 -maxcell 42 \
        R000466.head R000466.nx.hdf R7727_dark_062510.nx.hdf

which runs

    dtfind  R.head -cor flat.h5 -type 3 -seq 0 95 -pad 4 -window 20 20 -sigma 4 -min 5 -rect 2 12 253 1910
    dtindex Rfind.head Rfind.ref -reso 4 1 -known -threshold 100 -spacegroup N -fourier -errors 5 5 5 3 3 3 -nodeice
    dtrefine ... +CrysRot -wave 7.0 0.6 -rej 50 50 50 -go
    dtrefine ... +CrysRot +DetTrans1 +DetTrans2 -rej 1.5 1.5 1 -go
    dtrefine ... +CrysRot +DetAll -rej 1 1 1 -go
    dtintegrate Rrefine.head -flat -cor flat.h5 -type 204 -seq 0 95 -reso 1.0 20 -minpeakrad 2.0 \
                -h 100 -100 -k 100 -100 -l 100 -100 -rect 0 18 255 1910 -wave 7.0 0.6 -pad 4 -window 20 20
    dtref2afout -scale 1000 Ri.ref
    dtref2afout -scale auto -theta -keepneg -out Ri_aimless.afout Ri.ref
    dtexport -batch 466 -out Ri.mtz Rrefine.head Ri_aimless.afout

The last two steps (only if `dtexport` was built, which needs gemmi) write
`Ri.mtz`, an unmerged MTZ file for AIMLESS with one batch, numbered by the run
number, and the crystal orientation (see
[Export to MTZ](#export-to-mtz-for-aimless)). `Ri.afout` keeps the legacy
`negative.awk` format.

TOF-Laue options:

| Program | Option | Meaning |
|---|---|---|
| dtfind, dtintegrate | `-cor file -type mode` | Detector response from a V/Nb run. Modes: `none`, `pixel`, `segment`, `full`, `spectrum`. Legacy `3` = `spectrum` (divides by the incident spectrum; dtfind default `pixel`); legacy `204` = `full` (per-pixel response only; dtintegrate default). |
| dtfind | `-pad n` | TOF channels added on each side of a merged spot for its 3D centroid. |
| dtindex, dtrefine | `-threshold I` | Minimum spot intensity. |
| dtrefine, dtintegrate, dtpredict | `-wave λ1 λ2` | Wavelength range. |
| dtrefine | `-rej dx dy dchan` | The third limit is in TOF channels. |
| dtrefine | `+TOFDelay` | Refine the TOF delay (`Detector1_ROTATION_DELAY`). L1 is not refined because it is degenerate with the cell. |
| dtintegrate | `-h/-k/-l max min`, `-rect`, `-flat` | Index limits, detector rectangle, accepted for compatibility. |
| dtintegrate | `-tofint sum\|shoebox` | `sum` (default): integration by summation, see below. |
| dtintegrate | `-tofradius r [in out]` | Peak radius in pixels (default 3.5) and background annulus (default r+2.5 to r+6). |

Header keywords added to the legacy set: `TOF_WAVELENGTH_RANGE`,
`TOF_WIDTH_MODEL a b c` (predicted width in channels = a + bλ + c/sinθ;
default 2 0 0.5, as the original program used), `Detector1_SPATIAL_PARALLAX`,
`DTREFINE_RMS_CHAN`.

### Merging a data set

    scripts/pcs_merge.sh -head R000466refine.head -o merged -wave "0.6 6.3" R*/R*i.ref

The steps:
- `dtlauenorm` fits I = k_run · N(λ) · J_hkl. N is piecewise linear in ln λ with a smoothness restraint; outliers are rejected. It writes the normalized, unmerged list, with one batch per run.
- `dtscaleaverage` merges that list and reports statistics.
- `tools/pcs/ref2mtz.py` writes an unmerged MTZ with a LAMBDA column for AIMLESS.

For the original route, convert the lists with `dtref2afout -scale 1000` and
run LAUENORM. `dtref2afout` writes the same files as the legacy
`negative.awk`, with columns found by name. With the sigmas of this program,
use an I/σ limit of about 2 for LAUENORM's scale factors instead of the 7 in
the legacy `laue.dat`.

To scale with AIMLESS (including its Laue version), convert afout files to an
unmerged MTZ with `dtexport` (see [Export to MTZ](#export-to-mtz-for-aimless)).

### Model and conventions

- **Wavelength.** λ = (t(c) + delay) / (252.7784 µs m⁻¹ Å⁻¹ · (L1 + L2)).
  - t(c) is interpolated linearly between the TOF bin edges.
  - L2 is the flight path to the detector cylinder (radius = detector distance).
  - The observed wavelengths of the original dtfind are reproduced to < 0.0005 Å (`tCTOFLaue`).
- **Prediction.** For each hkl with d > d_min, λ = 2(X·S0)/|X|². If λ is within range, the spot position comes from the standard detector code and the channel from the TOF relation.
- **Lorentz factor.** `fCalc_lorentz` = 1/sin²θ, and `fIntensity` is divided by it, as in the original program. The remaining λ⁴, incident spectrum and detector efficiency are left to dtlauenorm or LAUENORM. LAUENORM applies no Lorentz correction to afout input.
- **Integration.** dtintegrate predicts every reflection of the run at once and integrates each from the in-memory volume:
  - the peak is a disk of radius r over the predicted TOF window, with fractional end channels;
  - the background is the mean of an annulus in each channel;
  - sigmas come from counting statistics.

  This integrates noticeably better than the rotation-image shoebox method (`-tofint shoebox`) and than the original program (merging statistics below).
- **Indexing and refinement.** These use reciprocal vectors X = (S + S0)/λ_obs. The radial (wavelength) part of the residual is down-weighted, and the residuals are reported in mm and in TOF channels.

### Differences from the original program

- **Intensity scale.** Integrated intensities are about 1100× smaller than the original program's (whose `-type 204` scaling is unknown), and dtfind intensities about 8× smaller. The legacy recipe thresholds (`dtindex -threshold 1000–4000`, `dtrefine -threshold 100`) therefore become about 100 and 10. LAUENORM rescales its input, so the scale doesn't matter there.
- **Sigmas.** The original program's sigmas are about 3.5× too small: its I/σ is inflated. The sigmas here are close to counting statistics. A legacy recipe's `-sigma 4` in dtfind therefore corresponds to about 1σ here.
- **Predicted channels.** These follow the cylinder flight path. The original predictions were up to about 0.4 channel later at the ends of the detector arc. Its observed wavelengths, however, used the cylinder path.
- **A legacy data-set error.** In the Crambin set, `R000588refine.head` points to `R000589.nx.hdf`, so the legacy "R000588" list is a second integration of R000589.

### Validation against the original outputs

| Check | Result |
|---|---|
| Stitching (2008, 2011, 2013 runs) | Median offset ≤ 0.06 px from legacy dtfind spots. |
| `tCTOFLaue` (λ ↔ channel, bin edges, legacy dtfind λ) | Pass on 5 runs (rms ≤ 0.0005 Å). |
| Prediction from a legacy header | rms 0.08 px, 0.0001 Å; all legacy reflections predicted. |
| dtfind (CA_Thiophene R000931) | 87% of all legacy spots and 91% of those with a valid width recovered. Centroids agree to 0.2 px and 0.19 channel. |
| dtindex + dtrefine, own pipeline (Crambin R000466) | Orientation within 0.06° of the legacy refinement; 0.54 mm, 0.39 channel. |
| dtintegrate vs legacy intensities (legacy I/σ > 3) | CC 0.966 mean over 47 Crambin runs; 0.962 CA R000931. GTA R7812: 0.82, which is 0.96 where both predict the same channel. |
| dtref2afout vs legacy `.afout` | 176 of 177 byte-identical; 1 legacy file was hand-edited later. |
| LAUENORM on identical reflections (Crambin, 47 runs) | RFACT1 0.238 here vs 0.249 for the legacy lists. |
| SCALA, legacy keywords (30–1.1 Å) | Rpim 0.151 vs 0.146; ⟨I/σ⟩ 3.6 vs 3.5; multiplicity 2.9 vs 2.8. |
| dtlauenorm on identical reflections | R 0.173 here vs 0.217 for the legacy lists. |
| X-ray regression (insulin, ADSC) | dtfind identical; dtindex same cell and orientation. |

### Not done

- Refinement of L1 (only the TOF delay is refinable).
- A physical depth model for the parallax. The legacy correction (vertical scale 1.006472 + 0.001257λ) is used.
- Automatic peak radius from the strong spots, and exclusion of neighbouring spots from the background annulus in the summation integration.
- Absorption and extinction corrections.

## Panel instruments (SNS MANDI, NXlauetof)

d\*TREK also reads TOF-Laue histograms of instruments with many flat
detector panels, in the NXlauetof layout of the ESS NMX reduced data, and
histograms SNS event data into that layout. Tested on SNS MANDI (CuZnSOD,
IPTS-37773, 10 runs, 40 panels of 256 × 256 pixels).

### Data files

- **NXlauetof** (`/entry/definition` = `NXlauetof`). One `NXdetector` per
  panel under `/entry/instrument`, with `data` (y, x, tof and an `axes`
  attribute), `origin` (centre of the first pixel), `fast_axis` and
  `slow_axis`, `x_pixel_size`, `y_pixel_size` and the TOF bin edges
  (`original_time_edges`, or `time_of_flight` centres). L1 is
  `/entry/instrument/source/distance`; the goniometer is the
  `depends_on` chain of `/entry/sample/transformations` (default: Mantid's
  omega, chi, phi). All panels must share the TOF bins.

  The panels are stacked along px1 into one mosaic image per TOF channel,
  separated by 8 masked rows; px0 is the panel's x index. The header gets
  `Detector1_SPATIAL_DISTORTION_TYPE=Panels_spatial`,
  `Detector1_PANELS_NUMBER` and one `Detector1_PANEL_k = bank start nx ny
  ox oy oz ux uy uz vx vy vz du dv` per panel (mm, d\*TREK laboratory frame).
  The detector goniometer is the identity, so the "mm" coordinates of a pixel
  are its laboratory position. `TOF_WIDTH_MODEL` defaults to
  `1 0 0.0075/r` for logarithmic bins of ratio r (1 0 1.5 for 0.5% bins).

- **SNS event NeXus** (`MANDI_*.nxs.h5`). `dteventbin` histograms the events
  in C++ (MANDI: 83 million events in 13 s, 70 MB output):

      dteventbin MANDI_13378.nxs.h5 -detcal calibration.DetCal \
                 -goniocal calibration_goniometer.xml [-out file] \
                 [-logbin 0.005 | -linbin dt] [-tof tmin tmax] [-bin n] \
                 [-gonio omega chi phi] [-banks 1,2,...]

  Panel geometry comes from the Mantid instrument definition in the file,
  or from an ISAW/Mantid `.DetCal` calibration (L1, T0 and panel positions)
  if one is given. The goniometer angles are time averages of the DAS logs.
  The goniometer calibration (chi offset, tilt) is recorded in the output but
  not applied; per-run refinement absorbs it.

- **`dtisawub UB.mat in.head out.head`** puts a Mantid/ISAW UB matrix into a
  header (for checks against Mantid).

### Processing and merging

    scripts/lauetof_process_run.sh -o r13378 -cell "65.73 65.73 151.88 90 90 120" -sg 179 \
        -detcal 2026B/calibration.DetCal -goniocal 2026B/calibration_goniometer.xml \
        MANDI_13378.nxs.h5

runs `dteventbin` (for event files), `dtfind` over all channels and the whole
mosaic, `dtindex -known -fourier`, three `dtrefine` stages and `dtintegrate
-tofradius 5`, writes the LAUENORM input with `dtref2afout -scale auto
-theta -keepneg`, and, if `dtexport` was built, `<run>i.mtz` for AIMLESS (one
batch, numbered by the run number unless `-batch n` is given, with the crystal
orientation). About 40 s per MANDI run. The last refinement stage refines an
unconstrained cell (`+CrysAll -CrysConstraints`), as Mantid refines an
unconstrained UB: it absorbs small systematic errors of the panel geometry and
places the predictions on the outer panels, which have few found spots,
2–5 px better.

To merge with LAUENORM, using for example the control file that Mantid's
`mandi_createmtz.py` writes:

    scripts/lauetof_lauenorm.sh -o merged -setup ~/Daresbury_laue/laue-setup.sh \
        -out merged.mtz lnorm_all_merged.dat r*/MANDI_*i.ref

or natively, `scripts/pcs_merge.sh -head r13378/r2.head ...` (`r2.head` has
the space-group-constrained cell), or with AIMLESS after `dtexport` (next
section).

In dtrefine, the mm residuals and the `-rej` limits of a panel detector are
along the axes of the panel of the predicted spot.

### Export to MTZ for AIMLESS

`dtexport` writes afout files (from `dtref2afout`, or Mantid's `laueNorm*`
files) as one unmerged MTZ file for AIMLESS, including its Laue version:

    dtexport -head r13378/r2.head -out mandi.mtz r*/refine.head r*/MANDI_*i.afout
    dtexport -cell 65.73 65.73 151.88 90 90 120 -spacegroup P6522 -out mantid.mtz laueNorm3*

- Each file becomes one batch, numbered from `-batch n` (default 1) in the
  order given. The batch headers carry the cell, the mean wavelength and
  `LDTYPE = 3` (Laue).
- The columns are `H K L M/ISYM BATCH I SIGI LAMBDA ROT`. The indices are
  reduced to the asymmetric unit with M/ISYM, and the file is sorted, as
  POINTLESS writes unmerged data. ROT is 0, so that AIMLESS does not take the
  batch number as the rotation angle (`-norot` omits it).
- afout files hold no cell or space group. They come from a d\*TREK header
  (`-head`) or from `-cell` and `-spacegroup`. The cell is made to obey the
  space group, with a warning if that changes it; use a header with the
  constrained cell (`r2.head` of `lauetof_process_run.sh`).
- Other options: `-wave lmin lmax`, `-scale f`, `-project`, `-crystal`,
  `-dataset`, `-title`. Reflections with σ ≤ 0 are skipped.
- AIMLESS needs the negative intensities, so write the afout files with
  `dtref2afout -keepneg`, as both scripts do. LAUENORM rejects zero and
  negative intensities alike, so its results do not change. `dtexport` warns
  if a file looks as if its negative intensities were set to 0.
- **Crystal orientation.** d\*TREK headers among the input files (the refined
  header of each run, one per afout file and in the same order) put the
  orientation of each run into its batch header, so that AIMLESS can apply its
  secondary-beam and absorption corrections. The header frame is the d\*TREK
  laboratory frame:
  - `UMAT` is the rotation U with U·B = G·(UB)_d\*TREK, where B is AIMLESS's
    (Busing & Levy) B matrix of the MTZ cell and G is the crystal goniometer
    matrix of the run. U is made orthogonal, which absorbs a refined cell that
    does not obey the space group (`refine.head` of `lauetof_process_run.sh`).
  - `SOURCE` is the d\*TREK source vector S0, anti-parallel to the beam.
  - `E1` is the outermost goniometer axis, with one axis at datum 0.
  - The rotation range is a nominal −0.5 to 0.5°, because AIMLESS ignores the
    orientation of a batch without one; ROT is 0.

  `dtexport` checks the wavelengths each header predicts for the reflections
  of its afout file. A header paired with the wrong file (MANDI: 24 % rms) is
  an error and nothing is written. On MANDI, the diffracted-beam directions
  that the MTZ headers give agree with the predicted pixel positions to a
  median 0.2° (maximum 0.45°).

  AIMLESS makes each oriented batch a run of its own (a stationary exposure),
  so `SCALES BATCH` fails ("no useful batch information"). Use
  `SCALES CONSTANT`, one scale per run, which gives the same result, and add
  `SECONDARY` or `ABSORPTION` as wanted.

`dtexport` uses the [gemmi](https://gemmi.readthedocs.io) C++ library
(`GEMMIHOME` at the top of the master `Makefile`, default
`/opt/homebrew/opt/gemmi`; `brew install gemmi`). `make check` lists the
libraries the build needs and fails if one is missing; `makeall.sh` runs it
before it cleans and rebuilds everything.

On the nine MANDI runs, the Laue version of AIMLESS (Chebyshev λ
normalization, 2–4 Å, batch scales, no orientation) gives the following:

| Statistic | This build | Mantid's laueNorm files |
|---|---|---|
| CC½ overall / outer shell (2.61–2.50 Å) | 0.978 / 0.461 | 0.971 / 0.478 |
| Rpim, all I+ and I− | 0.147 | 0.076 |
| Completeness | 99.3 % | 99.2 % |
| Observations / unique | 107,999 / 7,194 | 76,081 / 7,196 |

The two merged data sets correlate at 0.970 (0.61 at 2.50–2.66 Å). Mantid's
files hold 29 % fewer observations and no zero or negative intensities,
because Mantid drops weak reflections before writing them; this build writes
every predicted reflection. Keeping the weak ones raises Rpim, but CC½ is about
the same.

With the orientations (same keywords, `SCALES CONSTANT BFACTOR OFF`):

| Scale model | Rmerge | Rpim (inner shell) | ⟨I/σ⟩ | CC½ |
|---|---|---|---|---|
| No orientation, `SCALES BATCH` | 0.575 | 0.147 (0.051) | 5.3 | 0.978 |
| `CONSTANT`, `SECONDARY 0` | 0.575 | 0.147 (0.051) | 5.3 | 0.977 |
| `CONSTANT`, `SECONDARY 4` | 0.561 | 0.143 (0.046) | 5.6 | 0.980 |
| `CONSTANT`, `ABSORPTION 4` | 0.570 | 0.145 (0.053) | 5.4 | 0.979 |

### Validation against Mantid (MANDI CuZnSOD)

| Check | Result |
|---|---|
| `dteventbin` vs an independent histogram | Identical counts (banks 1 and 52); panel origins match the DetCal to < 1 µm. |
| Geometry vs Mantid | L2 rms 3 µm, scattering angles to 1e-5 rad, λ to 1e-4 Å. |
| Prediction with Mantid's UB (`dtisawub`) | 13017 predicted vs 13110 Mantid peaks; positions rms 0.29 px, TOF 0.001%. |
| Integration with Mantid's UB | CC 0.958 vs Mantid profile-fitted intensities (I/σ > 3, 3053 reflections). |
| Own pipeline, per run (`lauetof_process_run.sh`) | CC 0.900–0.966 vs Mantid (I/σ > 3; 2257–3952 reflections per run, matched by position). Refinement 0.77–0.89 mm, 0.61–0.82 channel. See below. |
| Merged, LAUENORM with Mantid's control file (9 runs, as Mantid) | CC(F²) 0.961 vs Mantid's merged MTZ (6573 reflections; 0.98 to 3.7 Å, 0.32 at 2.5–2.7 Å). RFACT1 0.172 vs Mantid's 0.214, on fewer accepted measurements (39 876 vs 62 267). |
| Merged natively (`pcs_merge.sh`, 9 runs) | CC 0.943 vs Mantid (7191 reflections; 0.52 at 2.5–2.7 Å). |

Per run (CC with Mantid's profile-fitted intensities for Mantid I/σ > 3):

| Run | Spots found | Refinement mm / channel | Integrated | CC (n) |
|---|---|---|---|---|
| 13378 | 405 | 0.89 / 0.82 | 12131 | 0.964 (2895) |
| 13379 | 427 | 0.82 / 0.80 | 12093 | 0.941 (3206) |
| 13380 | 486 | 0.81 / 0.75 | 12042 | 0.945 (3276) |
| 13382 | 429 | 0.88 / 0.74 | 12096 | 0.947 (2257) |
| 13383 | 512 | 0.84 / 0.69 | 12082 | 0.930 (3388) |
| 13384 | 622 | 0.82 / 0.61 | 12082 | 0.943 (3664) |
| 13385 | 710 | 0.80 / 0.74 | 12022 | 0.900 (3933) |
| 13386 | 683 | 0.77 / 0.64 | 12017 | 0.966 (3952) |
| 13387 | 657 | 0.81 / 0.72 | 12098 | 0.963 (3931) |
| 13388 | 636 | 0.83 / 0.71 | 12134 | 0.962 (3632) |

Our predictions and Mantid's agree to about 1 px (median) on most runs. In
run 13382 they differ by 2–4 px on the outer panels, which have almost no
found spots. There our predictions fit the found spots better (1.5 px rms)
than Mantid's peak positions do, so only 7965 of its reflections fall within
the 2.5 px matching window.

### Not done (panel instruments)

- The goniometer calibration is recorded, not applied. Each run's orientation is refined on its own, which absorbs it.
- No per-panel geometry refinement. The unconstrained cell absorbs small systematic errors.
- Weak reflections: summation integration gives a lower I/σ than Mantid's profile fitting, so LAUENORM rejects more of them (I < 1.5σ).
- Reflections on panel edges and in the gaps between panels are rejected.
