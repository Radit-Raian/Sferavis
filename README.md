# spherex_mover

Config-driven pipeline that finds **unlabeled moving objects** (candidate asteroids / comets) in
[SPHEREx](https://spherex.caltech.edu/) images. It downloads cutouts from IRSA around a sky position,
detects sources, removes known stars (Gaia DR3 / SIMBAD) and stationary sources, links the rest into
linear-motion tracklets, checks them against known Solar System bodies (SkyBoT) and writes a ranked
candidate list with diagnostic plots.

## Install

Requires **Python 3.10+** and an internet connection (IRSA, Gaia/VizieR, SkyBoT are queried live).

```bash
git clone <your-repo-url> spherex_mover
cd spherex_mover

python -m venv .venv
source .venv/bin/activate            # Windows (cmd):  .venv\Scripts\activate
                                     # Windows (PowerShell):  .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If a future library release breaks something, install the exact tested versions instead:
`python -m pip install -r requirements-tested.txt`.

## Run

```bash
python run_pipeline.py                              # asks for RA, Dec and search radius (deg)
python run_pipeline.py --ra 180 --dec 0 --radius 0.1
python run_pipeline.py --ra 20h23m57.8s --dec +38d30m27s --radius 0.1
python run_pipeline.py --no-prompt                  # reuse the target already in config.json
```

Or run the steps one at a time (each needs the previous one to have finished):
`python step1_config.py`, `python step2_download_obs.py`, ... `python step7_check_results.py`.
The scripts can be started from any folder; every path is resolved relative to this repo.

## Configuration

* `config.default.json` is the tracked template. On first run it is copied to **`config.json`**
  (git-ignored), which the pipeline rewrites with your target, output paths and progress flags.
  Edit `config.json` to tune thresholds (`detection`, `filtering`, `tracklets`, ...). To reset,
  delete `config.json`.
* All paths are stored relative to the repo with `/` separators, so a config works on Windows,
  macOS and Linux. Nothing user-specific is hard-coded.
* Optional environment variables:

| Variable | Meaning | Default |
|---|---|---|
| `SPHEREX_OUTPUT_DIR` | Where results/FITS cutouts go (can be absolute, e.g. a big disk) | `output/` inside the repo |
| `SPHEREX_CONFIG` | Use a different config file | `config.json` inside the repo |
| `IRSA_BASE_URL` | IRSA host used for cutout URLs | `https://irsa.ipac.caltech.edu/` |

## Steps

| Step | File | What it does |
|---|---|---|
| 1 | `step1_config.py` | validates the position, creates `output/<name>/` folders |
| 2 | `step2_download_obs.py` | IRSA TAP query + validated, resumable cutout download |
| 3 | `step3_detect_sources.py` | background, detection + deblending (photutils), pixel -> RA/Dec via WCS |
| 4 | `step4_survey_match.py` | Gaia DR3 (proper-motion corrected per epoch) and optional SIMBAD match |
| 5 | `step5_list_unlabeled_obj.py` | quality cuts, per-epoch Gaia astrometric calibration, static / moving classification, known-body recovery check |
| 6 | `step6_link_tracklets.py` | linear-motion tracklet linking within each visit |
| 7 | `step7_check_results.py` | SkyBoT known-object check, `candidates.csv`, ranked triage table, plots |
| - | `cutout_check.py` | optional diagnostic: shape / valid-pixel summary of the downloaded cutouts |

## Output

Everything for one target lives in `output/RA_<ra>_DEC_<dec>/` (images, `catalogs/`, `tracklets/`, `plots/`).
The main results are `tracklets/candidates.csv`, `tracklets/candidates_triage.csv` and `plots/candidate_*.png`.
`output/` is git-ignored: downloaded FITS cutouts can be large.

## Notes and troubleshooting

* SPHEREx images of one field are at different wavelengths, so brightness changes between epochs.
  Tracklets are linear fits, so keep `tracklets.visit_gap_days` small.
* Step 2 reuses `catalogs/cutout_urls.csv` and already-downloaded FITS files, so you can safely re-run it
  after an interruption. Delete that CSV if you change the cutout size or time range for the same target.
* `Missing Python packages` -> activate your virtual environment and run the install command above.
* Gaia download fails -> retry later, supply your own CSV (`ra, dec[, pmra, pmdec]`) in
  `config.json -> gaia.catalog_file`, or set `gaia.enabled` to `false`.
* On Windows, if a console shows garbled characters, run `set PYTHONUTF8=1` first.
* Candidates still need follow-up before any discovery claim.