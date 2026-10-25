# AGENTS.md — rules for any agent working in this repository

Read this file, the module's README and data/results/INDEX.md before any
task. State in your first message which existing result your task builds
on, or that it builds on none.

## Layout
- characterization/   pure Python: acceptance, solid angle, tolerances,
                      forward model (F), figures of merit. Reads meiga
                      output JSON and DetectorList.xml; defines no geometry
                      of its own beyond DetectorGeometry parameters.
- inversion/          FUENTES: 3D Bayesian inversion (Python, notebooks).
- thesis_orchestrator/ stages and contracts that call characterization and
                      inversion as black boxes and build the thesis report.
- sim/hodoscope-g4/   standalone Geant4 model of one module with optical
                      transport. Its geometry is the file
                      src/DetectorConstruction.cc.
- external/meiga/     submodule, branch hodoscope-na64mini: the two-module
                      telescope inside the Meiga framework. Its geometry is
                      src/G4Models/Hodoscope.cc + DetectorList.xml.
- config/hardware/    channel map and SiPM bias voltages (xlsx).
- docs/               hardware notes, KNOWN_ISSUES.md.
- data/raw/            measured raw inputs, outside git.
- data/ -> /mnt/d/muography-data, outside git.

## Where things go
- Source, Geant4 macros, LaTeX sources, configs: inside the owning module.
- Measured data (detector readouts, topography, test-beam reconstructions) lives in data/raw/<instrument>_<year>/ with MANIFEST.sha256; it is never regenerable and never deleted.
- Every run writes to data/results/<module>/<YYYY-MM-DD>_<tag>/ and must
  contain manifest.json: command, git commit of the monorepo and of the
  submodules, inputs with sha256, date, physical parameters (D, d, nbars,
  offset, pitch, reflector, R425, PDE, n_events, ...). Never write run
  outputs into the module tree.
- Figures used in a talk, paper or thesis chapter live next to the .tex
  that uses them; their manifest names the run that produced them.
- Exploratory work goes to data/scratch/<YYYY-MM-DD>_<tag>/. Anything
  there older than 30 days may be deleted without asking.
- Hardware documentation goes to docs/hardware/ or config/hardware/.
- After a run, add or update its row in data/results/INDEX.md.
- The Overleaf Muographers proceeding is at
  `/mnt/c/Users/renem/Dropbox/Aplicaciones/Overleaf/Muographers - Proceeding/`.
  Set `MUOGRAPHY_OVERLEAF_DIR` to override this default for figure syncs.

## Never
- Commit files > 5 MB, *.root, *.npz, build/, compiled PDFs, photos.
- Create a new top-level directory, a new studies/<name>/ outside
  sim/hodoscope-g4/studies/, or a second copy of a geometry definition.
  Geometry lives in exactly one place per simulation.
- Leave a result without a manifest; a result without a manifest is
  treated as scratch.
- Rename or move a file another module reads (see DEPENDENCIES.md)
  without updating the reader in the same commit.
- Use absolute paths; resolve from MUOGRAPHY_ROOT / MUOGRAPHY_DATA or
  from the file's own location.
- Create directory names in upper case or two paths that differ only by
  case (the repo is cloned on NTFS and macOS).
- Change a physical parameter default (D, pitch, offset, R425, PDE, ...)
  without a line in KNOWN_ISSUES.md or the module README saying why.

## Naming
- Directories and files: lower case, underscores, no spaces, no accents.
- Run tags: short, descriptive, no dates inside the tag (the directory
  carries the date).
- Variables keep the canonical names: D, d, nbars, offset, pitch, nph,
  npe, p_eff, R425, sigma_x, sigma_y, sigma_theta.
