# Neutron TOF-Laue data from the LANSCE PCS

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
unchanged: dtfind and dtindex give identical results on X-ray images.

## Data files

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

## Processing one run

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

## Merging a data set

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

## Model and conventions

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

## Differences from the original program

- **Intensity scale.** Integrated intensities are about 1100× smaller than the original program's (whose `-type 204` scaling is unknown), and dtfind intensities about 8× smaller. The legacy recipe thresholds (`dtindex -threshold 1000–4000`, `dtrefine -threshold 100`) therefore become about 100 and 10. LAUENORM rescales its input, so the scale doesn't matter there.
- **Sigmas.** The original program's sigmas are about 3.5× too small: its I/σ is inflated. The sigmas here are close to counting statistics. A legacy recipe's `-sigma 4` in dtfind therefore corresponds to about 1σ here.
- **Predicted channels.** These follow the cylinder flight path. The original predictions were up to about 0.4 channel later at the ends of the detector arc. Its observed wavelengths, however, used the cylinder path.
- **A legacy data-set error.** In the Crambin set, `R000588refine.head` points to `R000589.nx.hdf`, so the legacy "R000588" list is a second integration of R000589.

## Validation against the original outputs

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

## Not done

- Refinement of L1 (only the TOF delay is refinable).
- A physical depth model for the parallax. The legacy correction (vertical scale 1.006472 + 0.001257λ) is used.
- Automatic peak radius from the strong spots, and exclusion of neighbouring spots from the background annulus in the summation integration.
- Absorption and extinction corrections.
