# Run-2 KD plots

Confirmed the additional-jet criteria in active ANATree/src/TreeLoop.cc: resolved lines 604/626 and merged lines 1133/1159. Two-jet VBF is calculated for >=2 additional jets; one-jet VBF is calculated only for exactly one. Reset values at lines 3382–3385 are negative sentinels.

Kept the preceding Run-2 kinematic selections, normalization, blinding and full ParticleNet-corrected DY weights. All tags, channels and eras combined separately per topology. Read 336 inputs once for three variables; never opened SR Data.

Valid values are finite and in [0,1]. Full CR Data counts: resolved 408775 for kd, 68123 for kd_jjvbf, 125427 for kd_jvbf; merged 8780, 2734, 3157 respectively. No CR Data event has both VBF discriminants valid, consistent with the exclusive calculation branches. No CR Data values exceed one. These counts are from CR only; no observed SR yields were read.

Saved per-input validity/yield audit and histogram cache. Independent validation checks exact valid-event weighted yields and sumw2, Data counting variance, all 336 inputs, SR Data absence and 12 individual plus 4 overview PNG/PDF pairs. Existing Z kinematic plots remain untouched.
