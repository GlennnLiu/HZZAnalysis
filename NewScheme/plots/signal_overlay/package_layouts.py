"""Validate and package revised plots without modifying previous plot folders."""
import json,zipfile
from pathlib import Path

base=Path('/cms/user/liugeliang/CMS/HZZAnalysis/Plot')
folders=['run2_CR_m2l2q_DataMC','run2_SR_m2l2q_signals_layout','run2_discriminants_PNet_signals_layout',
         'run2_PNet_scores_signals_layout','run2_PNet_scores_linear_signals_layout',
         'run2_PNet_fullrange_signals_layout','run2_PNet_fullrange_linear_signals_layout',
         'run2_discriminants_PNet_layout','run2_PNet_scores_layout','run2_PNet_scores_linear_layout',
         'run2_PNet_fullrange_layout','run2_PNet_fullrange_linear_layout','run2_z_kinematics_PNet_layout']
counts={}; readme='''# Revised Run-2 plot layout, 2026-10-13

Legend and category labels remain inside the plotting axes; CMS labels and
luminosity remain at the original plot border. Legend side and y limits are
chosen using measured text boxes and the actual yields below those boxes.
Y limits follow the actual plotted yields and uncertainties, without large fixed
headroom multipliers. Original plots, trees, cards and input histogram caches
are unchanged. Revised plots retain original filenames in separate folders.

No signals are shown in CR, including the revised MELA and ParticleNet sets.
SR remains blinded; no SR Data file is loaded. Signals in SR remain separate
1 pb benchmark hypotheses with dotted ggF/solid VBF and mass-matched colors.
Nominal ParticleNet corrections are retained. The full-range ZvsQCD plots
continue to use no ParticleNet SF, following their existing convention.

New CR mass plots use direct DY MC + TTplusWW + VZ, with observed Data and a
Data/MC panel. MC statistical bands are sqrt(sumw2); Data errors are Garwood
68.27% Poisson intervals. CR is the existing full hadronic sideband, including
105–135 and VR; only CR files are read, without adding separate VR files.
All four eras retain original luminosity normalizations, displayed as 138 fb^-1.

Mass histograms share event-filled bins across tags/channels: 50 GeV from
300–1500, 100 GeV to 3500 and an open >=3500 bin displayed as 3500–4000.
Plots are density-scaled to Events / 50 GeV; overflow yields are retained.
Individual ee/mumu/combined and inclusive-tag plots are supplied. Ratio points
above 2.5 are indicated by triangles; nonpositive MC denominator bins are omitted
only from the ratio, not removed from the stored histograms.

Scripts and validation: ../signal_overlay_Run2_20261013/revise_layouts.py,
CR_mass_histograms.npz, CR_mass_input_yields.csv and layout_validation.json.
'''
for name in folders:
    folder=base/name
    png=list(folder.glob('*.png')); pdf=list(folder.glob('*.pdf'))
    expected=(28,29) if 'm2l2q' in name else ((20,21) if 'kinematics' in name else ((16,17) if 'discriminants' in name else ((3,4) if 'fullrange' in name else (6,7))))
    assert (len(png),len(pdf))==expected,(name,len(png),len(pdf),expected)
    assert all(p.stat().st_size>1000 for p in png+pdf)
    (folder/'README.md').write_text(readme)
    counts[name]=dict(png=len(png),pdf=len(pdf))
cache=base/'signal_overlay_Run2_20261013'
audit=json.loads((cache/'layout_validation.json').read_text())
assert len(audit)>100 and all(x['legend_inside_plot'] and x['measured_clearance'] for x in audit)
(cache/'revised_outputs.json').write_text(json.dumps(counts,indent=2))
with zipfile.ZipFile(base/'Run2_revised_layouts_and_CR_mass.zip','w',zipfile.ZIP_DEFLATED) as z:
    for name in folders:
        for f in (base/name).iterdir():
            if f.is_file(): z.write(f,f.relative_to(base))
    for name in ('revise_layouts.py','package_layouts.py','layout_validation.json','revised_outputs.json',
                 'CR_mass_histograms.npz','CR_mass_input_yields.csv','LAYOUT_WORKLOG.md'):
        z.write(cache/name,(cache/name).relative_to(base))
print('PASS:',len(folders),'folders,',sum(v['png'] for v in counts.values()),'PNGs,',sum(v['pdf'] for v in counts.values()),'PDFs; measured in-frame legend clearance.')
