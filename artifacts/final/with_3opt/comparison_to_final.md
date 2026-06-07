# Reduced 3-Opt Impact

This report compares the final MMAS experiments without local search against the final MMAS experiments with reduced 3-opt enabled.

## Findings

- **Reduced 3-opt average-best-deviation deltas: Baseline +4.22 pp, MSA heuristic +4.05 pp, Cycle cover +3.71 pp, Cycle-cover MSA patching +3.55 pp.**
- **Reduced 3-opt success-rate deltas: Baseline +31.85 pp, MSA heuristic +29.85 pp, Cycle cover +30.77 pp, Cycle-cover MSA patching +32.92 pp.**
- **Deviation gain over baseline with and without 3-opt: MSA heuristic +0.18 -> +0.01 pp, Cycle cover +0.58 -> +0.07 pp, Cycle-cover MSA patching +0.79 -> +0.12 pp.**
- **Signal remaining after enabling 3-opt: MSA heuristic 7.26%, Cycle cover 12.85%, Cycle-cover MSA patching 15.00%.**
- **This supports treating reduced 3-opt as a strong local-search layer that partially hides the construction heuristic effect.**

## Overall Effect

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>
</thead>
<tbody>
<tr><td>Baseline</td><td align="right">4.78</td><td align="right">0.56</td><td align="right">6.15</td><td align="right">38.00</td></tr>
<tr><td>MSA heuristic</td><td align="right">4.60</td><td align="right">0.55</td><td align="right">7.69</td><td align="right">37.54</td></tr>
<tr><td>Cycle cover</td><td align="right">4.20</td><td align="right">0.48</td><td align="right">4.92</td><td align="right">35.69</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">3.99</td><td align="right">0.44</td><td align="right">5.00</td><td align="right">37.92</td></tr>
</tbody>
</table>

## Heuristic Signal

<table>
<thead>
<tr><th>Heuristic</th><th>Dev. gain vs baseline without 3-opt [pp]</th><th>Dev. gain vs baseline with 3-opt [pp]</th><th>Signal remaining [%]</th></tr>
</thead>
<tbody>
<tr><td>MSA heuristic</td><td align="right">+0.18</td><td align="right">+0.01</td><td align="right">7.26</td></tr>
<tr><td>Cycle cover</td><td align="right">+0.58</td><td align="right">+0.07</td><td align="right">12.85</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">+0.79</td><td align="right">+0.12</td><td align="right">15.00</td></tr>
</tbody>
</table>

Deviation gain vs baseline is `baseline average best deviation - heuristic average best deviation`, so positive values mean that the heuristic improved over the baseline. Signal remaining is the share of this deviation gain still visible after enabling reduced 3-opt.
