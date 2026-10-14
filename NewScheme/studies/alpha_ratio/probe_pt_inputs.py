import json
from pathlib import Path
import sys
import uproot
import yaml

root = Path('/cms/user/liugeliang/CMS/HZZAnalysis')
sys.path.insert(0, str(root/'FitTemplate'))
import make_reduced_trees as producer
config = yaml.safe_load((root/'cards/config_UL2018.yml').read_text())
print('SELECTIONS', json.dumps(config['cut']['SR']['isEE'], indent=2), flush=True)
for sample in config['Samples_lists']:
    path = producer.input_path(config, '2018', sample)
    with uproot.open(path) as f:
        t = f['passedEvents']
        print(json.dumps(dict(sample=sample, path=str(path), entries=t.num_entries,
          pt_fields={k:k in t for k in ['pt2jet','ptmerged','EventWeight','w_deepjet']},
          norm=None if sample=='Data' else float(f['sumWeights'].values().sum()))), flush=True)
path = root/'FitTemplate/RedTrees_old/DY_ee_resolved_untag_2018_CR.root'
with uproot.open(path) as f:
    print('REDUCED_FIELDS', f['tree'].keys(), flush=True)
