# multijet_skim

Multijet ⟨x_j⟩ data/MC bias study (Oct 2026), which led to the selection and recoil
smearing now in `multijet/analysis.cc` (see `multijet/README.md`).

- `slides/multijet_slides.pdf`: skim mismatch, Jet5/Jet8 statistics, JER template, cut scans.
- `recoil_configs/slides/recoil_configs.pdf`: recoil-smearing configurations 1-4 and hybrid.
- `recoil_configs/config4_newmatch_truth7/`: the adopted configuration (truth2 matching,
  truth jets >= 7 GeV). R=0.4 data/MC ⟨x_j⟩: 0.999, 0.995, 0.981, 0.988 for leading pT
  20-25, 25-30, 30-35, 35-50 GeV. `config4_newmatch_truth5/` (5 GeV truth jets) is a
  possible systematic: 1.006, 0.999, 0.990, 0.993.

ROOT outputs (`*.root`, ~1.9 GB) are not in git; they stay on the laptop checkout. The
adopted configuration's `multijet_analysis_pythia.root` is also in Alpine's `multijet/`.
