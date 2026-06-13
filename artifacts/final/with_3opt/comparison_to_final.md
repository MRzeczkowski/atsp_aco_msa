# Reduced 3-Opt Impact

This report compares the final MMAS experiments without local search against the final MMAS experiments with reduced 3-opt enabled.

## Findings

- **Reduced 3-opt average-best-deviation deltas: Baseline +4.06 pp, MSA heuristic +3.95 pp, Cycle cover +2.53 pp, Cycle-cover MSA patching +1.78 pp.**
- **Reduced 3-opt success-rate deltas: Baseline +30.74 pp, MSA heuristic +34.37 pp, Cycle cover +27.19 pp, Cycle-cover MSA patching +28.30 pp.**
- **Deviation gain over baseline with and without 3-opt: MSA heuristic +0.14 -> +0.02 pp, Cycle cover +1.81 -> +0.28 pp, Cycle-cover MSA patching +2.54 -> +0.26 pp.**
- **Signal remaining after enabling 3-opt: MSA heuristic 15.85%, Cycle cover 15.24%, Cycle-cover MSA patching 10.29%.**
- **This supports treating reduced 3-opt as a strong local-search layer that partially hides the construction heuristic effect.**

## Overall Effect

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>
</thead>
<tbody>
<tr><td>Baseline</td><td align="right">4.60</td><td align="right">0.54</td><td align="right">9.56</td><td align="right">40.30</td></tr>
<tr><td>MSA heuristic</td><td align="right">4.46</td><td align="right">0.52</td><td align="right">9.78</td><td align="right">44.15</td></tr>
<tr><td>Cycle cover</td><td align="right">2.79</td><td align="right">0.26</td><td align="right">7.41</td><td align="right">34.59</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">2.06</td><td align="right">0.28</td><td align="right">7.85</td><td align="right">36.15</td></tr>
</tbody>
</table>

## Heuristic Signal

<table>
<thead>
<tr><th>Heuristic</th><th>Dev. gain vs baseline without 3-opt [pp]</th><th>Dev. gain vs baseline with 3-opt [pp]</th><th>Signal remaining [%]</th></tr>
</thead>
<tbody>
<tr><td>MSA heuristic</td><td align="right">+0.14</td><td align="right">+0.02</td><td align="right">15.85</td></tr>
<tr><td>Cycle cover</td><td align="right">+1.81</td><td align="right">+0.28</td><td align="right">15.24</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">+2.54</td><td align="right">+0.26</td><td align="right">10.29</td></tr>
</tbody>
</table>

Deviation gain vs baseline is `baseline average best deviation - heuristic average best deviation`, so positive values mean that the heuristic improved over the baseline. Signal remaining is the share of this deviation gain still visible after enabling reduced 3-opt.
