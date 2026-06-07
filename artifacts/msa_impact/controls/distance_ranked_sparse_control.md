# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.40`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 6/8 instances.**
- **Mean average best deviation: MSA 1.80%, distance-ranked sparse 2.54%, delta -0.74 pp.**
- **Mean success rate: MSA 17.92%, distance-ranked sparse 15.42%, delta +2.50 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.125000.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">2.78</td><td align="right">-0.66</td><td align="right">3.33</td><td align="right">13.33</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">100.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.77</strong></td><td align="right">5.19</td><td align="right">-1.42</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.09</strong></td><td align="right">1.22</td><td align="right">-0.13</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv64</td><td align="right">0.99</td><td align="right"><strong>0.86</strong></td><td align="right">+0.13</td><td align="right">40.00</td><td align="right">10.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.48</strong></td><td align="right">5.42</td><td align="right">-2.94</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.24</strong></td><td align="right">3.76</td><td align="right">-0.52</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.74</strong></td><td align="right">1.09</td><td align="right">-0.35</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.80</strong></td><td align="right"><strong>2.54</strong></td><td align="right"><strong>-0.74</strong></td><td align="right"><strong>17.92</strong></td><td align="right"><strong>15.42</strong></td></tr>
</tbody>
</table>
