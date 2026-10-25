# Known Issues and Physical Parameter Discrepancies (`docs/KNOWN_ISSUES.md`)

This document records identified physical parameter conventions and geometric discrepancies between the simulation frameworks in the consolidated repository. **Neither geometry definition has been modified** in Phase 1b to preserve baseline reproducibility.

---

## 1. Sub-Plane Z-Spacing Discrepancy (3.0 mm vs 2.0 mm)

The internal Z-spacing between paired subplanes (subplane 0/A and subplane 1/B of the same orientation) differs between `sim/hodoscope-g4` and `external/meiga`:

* **`sim/hodoscope-g4` (3.0 mm centre-to-centre):**
  - **Location:** [`sim/hodoscope-g4/src/DetectorConstruction.cc:117-119`](file:///home/reriosto/muography-thesis/sim/hodoscope-g4/src/DetectorConstruction.cc#L117-L119) and [`sim/hodoscope-g4/src/DetectorConstruction.cc:370`](file:///home/reriosto/muography-thesis/sim/hodoscope-g4/src/DetectorConstruction.cc#L370).
  - **Definition:**
    ```cpp
    fBarThickness ( 1.0 * mm),
    fSepZ         ( 2.0 * mm),
    ...
    G4double dz_sup_inf = fBarThickness + fSepZ; // 1.0 mm + 2.0 mm = 3.0 mm
    ```
  - **Interpretation:** `fSepZ` is treated as a physical surface gap between the 1.0 mm thick scintillator bars, resulting in a 3.0 mm centre-to-centre separation.

* **`external/meiga` (2.0 mm centre-to-centre):**
  - **Location:** [`external/meiga/src/G4Models/Hodoscope.cc:171-174`](file:///home/reriosto/muography-thesis/external/meiga/src/G4Models/Hodoscope.cc#L171-L174).
  - **Definition:**
    ```cpp
    const double localZ =
        (axis == "X")
            ? ((subplane == 0) ? 1.5 * CLHEP::mm : 3.5 * CLHEP::mm)
            : ((subplane == 0) ? 8.5 * CLHEP::mm : 10.5 * CLHEP::mm);
    ```
  - **Interpretation:** Subplane centres are explicitly placed at `1.5 mm` and `3.5 mm` for X, and `8.5 mm` and `10.5 mm` for Y, giving a centre-to-centre separation of `3.5 - 1.5 = 2.0 mm` (surface gap of 1.0 mm).

---

## 2. X–Y Plane Separation Convention (Gap vs Centre-to-Centre)

The physical distance between the X-measuring plane and the Y-measuring plane within a single station uses different geometric conventions:

* **`sim/hodoscope-g4` (Surface-to-Surface Gap):**
  - **Location:** [`sim/hodoscope-g4/src/DetectorConstruction.cc:368,371`](file:///home/reriosto/muography-thesis/sim/hodoscope-g4/src/DetectorConstruction.cc#L368-L371).
  - **Definition:**
    ```cpp
    // D es gap entre superficies
    G4double dz_X_Y = fBarThickness + fD;
    ```
  - **Interpretation:** The parameter `D` (default 5.0 mm via UI `/hodoscope/det/setD`) defines the clear air gap between the inner surfaces of X-inf and Y-sup. The centre-to-centre distance is therefore `fBarThickness + fD = 1.0 + 5.0 = 6.0 mm`.

* **`external/meiga` (Centre-to-Centre Distance):**
  - **Location:** [`external/meiga/src/G4Models/Hodoscope.cc:171-174`](file:///home/reriosto/muography-thesis/external/meiga/src/G4Models/Hodoscope.cc#L171-L174).
  - **Definition:**
    ```cpp
    (axis == "X") ? ... 3.5 * CLHEP::mm ... : ... 8.5 * CLHEP::mm ...
    ```
  - **Interpretation:** The centre of X-inf is at `3.5 mm` and the centre of Y-sup is at `8.5 mm`. The centre-to-centre distance between X and Y is fixed at `8.5 - 3.5 = 5.0 mm`, which corresponds to an air gap of `4.0 mm`.

### Resolved 2026-10-24

The measured-CAD reference is the STEP model **“Abrazadera Lateral superior der pt1”**: its comb rows are at `z = 4–5 mm` and `z = 7–8 mm`, with 4.00 mm slot pitch and slot centres spanning `-15..15 mm`. René's statement of 2026-05-13 confirms that 2 mm and 5 mm denote clear air gaps. Therefore both implementations now use a 3.0 mm same-axis centre spacing (2 mm face gap) and a 6.0 mm X-block/Y-block centre spacing (5 mm inner-face gap).

`sim/hodoscope-g4` is centred on module `Z=0` with layer centres `+6,+3,-3,-6 mm`; Meiga uses module-local coordinates `1.5,4.5,10.5,13.5 mm` within `0..15 mm`. Both use union-centred transverse coordinates A=`-15..13 mm`, B=`-13..15 mm`, at 4 mm pitch. Meiga's inter-module `D=98.4 mm` remains the centre-to-centre separation of its two modules.

The older paragraphs above describe the pre-resolution code and naming convention; the CAD-backed values above are the resolved geometry.

---

## 3. Impact and Resolution Policy
- Downstream reconstruction and characterization modules (`characterization/` and `inversion/`) rely on `effective_station_separation_mm` (nominal 98.4 mm) and define station positions through `characterization/forward/`.
- At the end of Phase 1b, the two detector geometries still differed. The CAD-backed reconciliation below was applied later on the geometry fix branch.

### Resolved 2026-10-24

The CAD and René's 2026-05-13 clarification resolve both discrepancies described above. The root branch `geom/fix-2026-10-24` now aligns both simulation geometries with the stated 3.0 mm and 6.0 mm centre spacings and 2.0 mm and 5.0 mm clear air gaps. The 98.4 mm inter-module centre separation remains unchanged.

Geometry outputs carry a version marker: version 1 identifies runs dated on or before 2026-10-23, which used the earlier geometry; version 2 identifies outputs produced with the CAD-aligned geometry from 2026-10-24 onward. hodoscope-g4 stores this as ROOT `TNamed` `geometry_version`; Meiga NA64Mini stores it as the top-level JSON key `geometry_version`.

## 4. EcoMug normalisation

A paired one-million-primary check was run on 2026-10-24 using EcoMug's `GetEstimatedTime(1000000)` (default 10-million-point integration), 0.1–10,000 GeV/c momentum, and identical seeds 20,461,024–20,461,063. The eSky surface was 1.0 × 1.0 m at z=+0.30 m; eHalfSphere used R=0.25 m. For both, the saved EcoMug generation point was preserved and each downward track was extrapolated to the centered 0.20 × 0.20 m square at z=0.

| Surface | Crossings | Exposure | Flux through square |
|---|---:|---:|---:|
| eSky | 37,928 / 1,000,000 | 8,098.71 s | 117.08 ± 0.59 (counting) ± 1.90 (integration) Hz/m² |
| eHalfSphere | 176,164 / 1,000,000 | 35,930.41 s | 122.57 ± 0.27 (counting) ± 2.34 (integration) Hz/m² |

The eHalfSphere/eSky flux ratio is 1.047, near unity and inconsistent with the proposed factor-of-two surface-area correction. eSky itself is 9.2% below the declared 129 Hz/m² after combining counting and integration errors, so this test does **not** establish 129 Hz/m² within statistics and does not isolate a hemisphere-versus-sky exposure error. The finite 1 m eSky plane is 0.30 m above the test square, so oblique tracks crossing the square can extrapolate outside the generated plane. No third re-scaling of the v2b rate table is justified by this test.

The earlier v2b crossing check is not a valid exposure-factor measurement: Meiga's production `PrimaryGenerator.cc` invokes `ComputeInjectionPosition` after EcoMug `Generate()`, replacing the generator's position with a separately sampled XML injection position. EcoMug also uses SI meters (`EcoMug.h`, `EMUnits`), whereas Meiga's XML dimensions are Geant4 lengths in mm. The Round 4 test used a temporary build that selected eSky/eHalfSphere, passed EcoMug geometry in meters, and converted the original EcoMug position to Geant4 mm without resampling; no source or submodule edit was made in the repository.

**Exposure convention going forward:** use EcoMug's own `GetEstimatedTime(n)` for the configured generation surface, passing its dimensions in meters and recording whether the method is eSky, eCylinder, or eHalfSphere. Preserve the generated EcoMug position when transporting the event, converting units only at the Geant4 boundary. Do not apply the Round 3 empirical 0.469 multiplier or substitute πR² for the documented eHalfSphere curved area 2πR². Treat the 129 Hz/m² closure as unresolved until the finite-surface and Meiga position/unit path are validated in the production generator.
