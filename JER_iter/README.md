# JER_iter: JES <-> JER iteration study

Iterative in-situ JES / JER with the jet group. We send an in-situ JES file, they send back JER
smearing templates, and we redo the in-situ JES with each new JER (nominal systag only). It runs the
main analysis code unchanged, with no separate branch or worktree:

- `GAMMAJET_JER` (`ana::mcJerMode`) sets the MC jet pT:
  - unset or `tree`: the tree's smeared branches. This is the official analysis.
  - `none`: `jet_pt_calib`.
  - a template file: re-smeared per jet from `jet_smear_z` / `jet_truthmatch_pt`. This needs the
    v1jerz Photon trees.
- `GAMMAJET_JES_TABLE` (`ana::jesForSystag`) replaces the compiled ana.h p_a table with a text file.
  Each pass updates the file, so ana.h is never edited and nothing is rebuilt.
- `GAMMAJET_UNFOLD` points main's macros at a per-run work dir, `work/<tag>/`. Main's `hists/`,
  `insitu/` and `reweight/` outputs are never touched.

The in-situ scan is main's `grid_insitu.C`, i.e. the MC-folded scan.

## Layout

| Path | Contents |
|---|---|
| `run_iter.sh TAG JER [START_TABLE]` | One iteration run: passes of unfold -> purity -> reweight -> unfold -> scan until p_a moves < 0.0005 (max 4). Writes `work/TAG/insitu_jes_TAG.root` |
| `setup_work.sh TAG` | Makes `work/TAG/` (copies of main's hists/reweight, `trees` -> `../../trees`) |
| `closure_jerfile.sh` | Closure: template re-smear (file mode) vs the tree's smeared branches |
| `read_pa.C`, `write_jes_file.C`, `write_table.C`, `compare_files.C` | Helpers |
| `inputs/` | Files received: `bisectorJER_insituJES1.root` (iter 1 JER, md5 443ccf53); `jes_table_iter1_standard.txt` (all systags at the iter-1 p_a, for the on-hold standardization) |
| `templates/` | `jerband_smearing_templates.root`, the JER the v1/v1jerz trees were smeared with |
| `out/` | Files sent back: `insitu_jes_iter0_noJER.root`, `insitu_jes_iter1_bisector.root` |
| `trees/` (not in git) | v1jerz pythia Photon5/10/20 trees, plus links to main's `trees/` for the rest |
| `work/` (not in git) | Per-run work dirs and logs; `work/archive_worktree/` holds the old PPG18_jeriter runs |

## Running

Long jobs (about 8 min per pass): run them under Slurm or in the background, e.g.

    JER_iter/run_iter.sh iter2_<name> /path/to/JER_templates.root
    GAMMAJET_JER_HIST='h_quaddiff_bisector_r{R}' JER_iter/run_iter.sh iter1_bisector JER_iter/inputs/bisectorJER_insituJES1.root

`{R}` in `GAMMAJET_JER_HIST` becomes the radius (`03` for R=0.3). A radius with no TH1 uses the
nearest one, with a warning; iter 1 has no r02, so R=0.2 uses r03.

## History

| Iteration | MC JER | p_a R=0.2 ... 0.8 | Notes |
|---|---|---|---|
| iter0 | none (`jet_pt_calib`) | 0.9280 0.9100 0.9283 0.9212 0.9166 0.9130 0.9246 | Oct 27, old scan |
| iter1 | `bisectorJER_insituJES1.root` | 0.9136 0.9058 0.9212 0.9170 0.9077 0.9116 0.9203 | Oct 29, old scan, converged in 2 passes |

Both used the pre-Nov 3 scan in the PPG18_jeriter worktree, where ana.h was edited and rebuilt.
Reruns here use the MC-folded scan, so expect small shifts.

The standardization (every systag at the iter-1 JES, MC re-skimmed with the bisector JER) is on hold
until the full JER errors arrive. `inputs/jes_table_iter1_standard.txt` holds its JES table.
