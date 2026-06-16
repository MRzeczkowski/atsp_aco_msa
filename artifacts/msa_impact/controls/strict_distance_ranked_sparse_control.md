# Strict MSA Distance-ranked Sparse Control

This sanity check compares Strict MSA against a deterministic sparse mask built from the cheapest directed edges. The control preserves the matching edge-count structure and uses the best Strict MSA impact heuristic weight selected separately for each instance.

## Findings

- **Strict MSA had lower average best deviation than the strict distance-ranked sparse control in 7/7 instances.**
- **Mean average best deviation: Strict MSA 1.67%, strict distance-ranked sparse 2.78%, delta -1.11 pp.**
- **Mean success rate: Strict MSA 10.00%, strict distance-ranked sparse 2.86%, delta +7.14 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means Strict MSA had lower average best deviation than the strict distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>Strict MSA avg best dev. [%]</th><th>strict distance-ranked sparse avg best dev. [%]</th><th>Delta [pp]</th><th>Strict MSA success [%]</th><th>strict distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">3.47</td><td align="right">-1.35</td><td align="right">20.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.25</strong></td><td align="right">5.31</td><td align="right">-3.06</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.03</strong></td><td align="right">1.13</td><td align="right">-0.10</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.36</strong></td><td align="right">0.72</td><td align="right">-0.36</td><td align="right">50.00</td><td align="right">20.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.10</strong></td><td align="right">4.03</td><td align="right">-1.93</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.32</strong></td><td align="right">3.53</td><td align="right">-0.21</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.53</strong></td><td align="right">1.27</td><td align="right">-0.74</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.67</strong></td><td align="right"><strong>2.78</strong></td><td align="right"><strong>-1.11</strong></td><td align="right"><strong>10.00</strong></td><td align="right"><strong>2.86</strong></td></tr>
</tbody>
</table>
