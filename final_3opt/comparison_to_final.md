# Reduced 3-Opt Impact

This report compares the final MMAS experiments without local search against the final MMAS experiments with reduced 3-opt enabled.

## Findings

- **Reduced 3-opt average-best-deviation deltas: Baseline +3.96 pp, MSA heuristic +2.92 pp, Cycle cover +3.01 pp, Cycle-cover MSA patching +3.06 pp.**
- **Reduced 3-opt success-rate deltas: Baseline +62.80 pp, MSA heuristic +55.60 pp, Cycle cover +61.20 pp, Cycle-cover MSA patching +58.40 pp.**
- **Deviation gain over baseline with and without 3-opt: MSA heuristic +1.07 -> +0.03 pp, Cycle cover +0.92 -> -0.03 pp, Cycle-cover MSA patching +0.93 -> +0.03 pp.**
- **Signal remaining after enabling 3-opt: MSA heuristic 2.63%, Cycle cover 3.46%, Cycle-cover MSA patching 3.02%.**
- **This supports treating reduced 3-opt as a strong local-search layer that partially hides the construction heuristic effect.**

## Overall Effect

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>
</thead>
<tbody>
<tr><td>Baseline</td><td align="right">4.09</td><td align="right">0.13</td><td align="right">2.80</td><td align="right">65.60</td></tr>
<tr><td>MSA heuristic</td><td align="right">3.02</td><td align="right">0.10</td><td align="right">9.60</td><td align="right">65.20</td></tr>
<tr><td>Cycle cover</td><td align="right">3.17</td><td align="right">0.16</td><td align="right">2.00</td><td align="right">63.20</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">3.16</td><td align="right">0.10</td><td align="right">6.80</td><td align="right">65.20</td></tr>
</tbody>
</table>

## Heuristic Signal

<table>
<thead>
<tr><th>Heuristic</th><th>Dev. gain vs baseline without 3-opt [pp]</th><th>Dev. gain vs baseline with 3-opt [pp]</th><th>Signal remaining [%]</th></tr>
</thead>
<tbody>
<tr><td>MSA heuristic</td><td align="right">+1.07</td><td align="right">+0.03</td><td align="right">2.63</td></tr>
<tr><td>Cycle cover</td><td align="right">+0.92</td><td align="right">-0.03</td><td align="right">3.46</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">+0.93</td><td align="right">+0.03</td><td align="right">3.02</td></tr>
</tbody>
</table>

Deviation gain vs baseline is `baseline average best deviation - heuristic average best deviation`, so positive values mean that the heuristic improved over the baseline. Signal remaining is the share of this deviation gain still visible after enabling reduced 3-opt.
