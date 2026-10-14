# All-era reduced trees with DY ParticleNet weights

Production initiated 2026-10-09. This directory retains the producer, streaming
runner, unit tests and batch scripts. The remote reproducible snapshot is:

`/cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/jobs/reduced_20261009_pnet`

New outputs (existing reduced trees are untouched):

`/cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet`

## Scope and weights

Four eras: 2016preAPV, 2016postAPV, 2017, 2018. Both topologies, three tags,
separate ee/mumu channels and CR/SR/VR. Processes: DY, TTplusWW, VZ, Data.
There are 72 ROOT files per era/topology, 576 in total. TTplusWW includes
TTTo2L2Nu, WWTo2L2Nu and tZq; VZ includes WZTo2Q2L and ZZTo2Q2L.

- `weight`: normalized full signed nominal weight, including applicable ParticleNet.
- `weight_no_pnet`: the same weight excluding only ParticleNet. Lepton and any
  applicable existing DeepJet corrections remain included.
- `pnet_sf`: product of `pnet_sf_zvsqcd` and `pnet_sf_zvslight`.
- `weight = weight_no_pnet * pnet_sf` event by event.

DY merged events use the era/flavour/pT/score tables in the frozen `sf/` directory,
not the legacy constant 1.36. These are the existing analysis calibrations,
approved by the user for DY, not a new independent certification of their physics
applicability. Merged VBF uses ZvsQCD; merged non-VBF also uses ZvsLight.
Resolved ParticleNet factors are one. DY DeepJet remains unapplied, matching
the existing workflow. tZq ParticleNet remains one, matching the existing helper.
Existing resolved non-DY DeepJet weights are used only in non-VBF categories;
their legacy efficiency implementation has not been changed in this production.
Negative MC weights are preserved; missing SF coverage raises an error.

## Features and cuts

Sixteen float64 branches: `mass`, `m2l2q`, `pt2l`, `mass2l`, `pt2q`, `mass2q`,
`kd`, `kd_jjvbf`, `kd_jvbf`, `particleNetZvsQCD`, `particleNetZbbvslight`,
and the five weight fields above. `mass` is retained as a compatibility alias
equal to `m2l2q`. Hadronic quantities and KDs map to the topology's native
source branches. ParticleNet scores in resolved trees remain the original
event's AK8 quantities; they are not a resolved-jet calibration.

Use the era card selections with an explicit lepton-channel guard. The latter
prevents missing `isEE` clauses in some older cuts from duplicating muons into ee.
Other era-specific selection choices are retained, including existing MET cuts.
CR is the card's `(m2q < 70) | (m2q > 105)`, including both 105–135 and VR.
SR has strict 70–105 bounds; VR has strict 135–165 bounds. Three older cards'
duplicated SR-as-VR expressions were corrected with backups retained.
No new hadronic 40–180 acceptance or m2l2q cut is imposed. A held-out VR alpha
calculation must explicitly exclude VR events from these full-CR trees.

## Validation and status

Read only one combined skim file per sample, never combined and split files
together. Normalize backgrounds once with the normal bins of `sumWeights`.
Processing is chunked at 100,000 entries; outputs remain `.partial` until every
tree in a job has passed a complete reread. A `validation_{era}_{topology}.json`
and `{era}_{topology}.done` certify each completed job. Reports record input
metadata, selection/config/code/SF hashes, entries and signed-weight statistics.

All eight era/topology input/schema/region preflights passed (11 inputs each).
Five unit tests passed. A 2018 merged pilot over 3,000 input entries per sample
passed all checks, writing 72 trees. Full production completed on 2026-10-09:
all eight jobs finalized and the aggregate validation passed for 576 ROOT trees. Absence
of `.done` means that job is not finalized. No alpha files or plots are regenerated.

| Era | Resolved job | Merged job |
| --- | --- | --- |
| 2016preAPV | 81981959 | 81981960 |
| 2016postAPV | 81981961 | 81981962 |
| 2017 | 81981963 | 81981964 |
| 2018 | 81981965 | 81981958 |

Check the scheduler using the full command if `hep_q` is not on PATH:

```bash
/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/hep_q -u liugeliang
```

Logs are `{era}_{topology}.out` and `.err` in the remote job directory.
Once all eight jobs finalize, run the read-only aggregate validator:

```bash
cd /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/jobs/reduced_20261009_pnet
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python validate_production.py ../../RedTrees_20261009_pnet
```

Do not resubmit into a partially populated output directory. The runner refuses
to overwrite existing final or partial trees. Use a fresh output directory for
an intentional replacement run. Data trees are retained for all regions, but
SR observations are not plotted by this workflow.
