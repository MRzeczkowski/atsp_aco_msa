# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the best MSA-impact heuristic weight selected separately for each instance.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 7/7 instances.**
- **Mean average best deviation: MSA 1.71%, distance-ranked sparse 3.11%, delta -1.41 pp.**
- **Mean success rate: MSA 11.43%, distance-ranked sparse 1.43%, delta +10.00 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">3.47</td><td align="right">-1.35</td><td align="right">20.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.61</strong></td><td align="right">3.59</td><td align="right">-0.98</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.06</strong></td><td align="right">1.13</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.36</strong></td><td align="right">1.12</td><td align="right">-0.76</td><td align="right">60.00</td><td align="right">10.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.08</strong></td><td align="right">6.08</td><td align="right">-4.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.06</strong></td><td align="right">5.04</td><td align="right">-1.98</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.65</strong></td><td align="right">1.35</td><td align="right">-0.70</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.71</strong></td><td align="right"><strong>3.11</strong></td><td align="right"><strong>-1.41</strong></td><td align="right"><strong>11.43</strong></td><td align="right"><strong>1.43</strong></td></tr>
</tbody>
</table>
