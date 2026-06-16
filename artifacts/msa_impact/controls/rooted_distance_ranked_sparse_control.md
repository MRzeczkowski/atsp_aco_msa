# Rooted MSA Distance-ranked Sparse Control

This sanity check compares Rooted MSA against a deterministic sparse mask built from the cheapest directed edges. The control preserves the matching edge-count structure and uses the best Rooted MSA impact heuristic weight selected separately for each instance.

## Findings

- **Rooted MSA had lower average best deviation than the rooted distance-ranked sparse control in 7/7 instances.**
- **Mean average best deviation: Rooted MSA 1.65%, rooted distance-ranked sparse 2.64%, delta -0.99 pp.**
- **Mean success rate: Rooted MSA 7.14%, rooted distance-ranked sparse 1.43%, delta +5.71 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means Rooted MSA had lower average best deviation than the rooted distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>Rooted MSA avg best dev. [%]</th><th>rooted distance-ranked sparse avg best dev. [%]</th><th>Delta [pp]</th><th>Rooted MSA success [%]</th><th>rooted distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>0.76</strong></td><td align="right">3.02</td><td align="right">-2.26</td><td align="right">10.00</td><td align="right">10.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.28</strong></td><td align="right">3.82</td><td align="right">-1.54</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.05</strong></td><td align="right">1.20</td><td align="right">-0.15</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.45</strong></td><td align="right">1.24</td><td align="right">-0.79</td><td align="right">40.00</td><td align="right">0.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>3.67</strong></td><td align="right">3.89</td><td align="right">-0.22</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right"><strong>2.67</strong></td><td align="right">4.10</td><td align="right">-1.43</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.70</strong></td><td align="right">1.21</td><td align="right">-0.51</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.65</strong></td><td align="right"><strong>2.64</strong></td><td align="right"><strong>-0.99</strong></td><td align="right"><strong>7.14</strong></td><td align="right"><strong>1.43</strong></td></tr>
</tbody>
</table>
