# Reduced 3-Opt Impact

This report compares the final MMAS experiments without local search against the final MMAS experiments with reduced 3-opt enabled.

## Findings

- **Reduced 3-opt average-best-deviation deltas: Baseline +3.16 pp, MSA heuristic +2.77 pp, Cycle cover +2.53 pp, Cycle-cover MSA patching +2.43 pp.**
- **Reduced 3-opt success-rate deltas: Baseline +34.07 pp, MSA heuristic +33.36 pp, Cycle cover +35.79 pp, Cycle-cover MSA patching +38.00 pp.**
- **Deviation gain over baseline with and without 3-opt: MSA heuristic +0.41 -> +0.02 pp, Cycle cover +0.74 -> +0.11 pp, Cycle-cover MSA patching +0.82 -> +0.09 pp.**
- **Signal remaining after enabling 3-opt: MSA heuristic 5.49%, Cycle cover 15.18%, Cycle-cover MSA patching 11.20%.**
- **This supports treating reduced 3-opt as a strong local-search layer that partially hides the construction heuristic effect.**

## Overall Effect

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>
</thead>
<tbody>
<tr><td>Baseline</td><td align="right">3.63</td><td align="right">0.47</td><td align="right">3.93</td><td align="right">38.00</td></tr>
<tr><td>MSA heuristic</td><td align="right">3.22</td><td align="right">0.45</td><td align="right">5.64</td><td align="right">39.00</td></tr>
<tr><td>Cycle cover</td><td align="right">2.89</td><td align="right">0.36</td><td align="right">4.36</td><td align="right">40.14</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">2.81</td><td align="right">0.38</td><td align="right">4.29</td><td align="right">42.29</td></tr>
</tbody>
</table>

## Heuristic Signal

<table>
<thead>
<tr><th>Heuristic</th><th>Dev. gain vs baseline without 3-opt [pp]</th><th>Dev. gain vs baseline with 3-opt [pp]</th><th>Signal remaining [%]</th></tr>
</thead>
<tbody>
<tr><td>MSA heuristic</td><td align="right">+0.41</td><td align="right">+0.02</td><td align="right">5.49</td></tr>
<tr><td>Cycle cover</td><td align="right">+0.74</td><td align="right">+0.11</td><td align="right">15.18</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">+0.82</td><td align="right">+0.09</td><td align="right">11.20</td></tr>
</tbody>
</table>

Deviation gain vs baseline is `baseline average best deviation - heuristic average best deviation`, so positive values mean that the heuristic improved over the baseline. Signal remaining is the share of this deviation gain still visible after enabling reduced 3-opt.
