# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the best MSA-impact heuristic weight selected separately for each instance.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 7/7 instances.**
- **Mean average best deviation: MSA 1.85%, distance-ranked sparse 2.58%, delta -0.74 pp.**
- **Mean success rate: MSA 5.71%, distance-ranked sparse 0.00%, delta +5.71 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>0.70</strong></td><td align="right">2.14</td><td align="right">-1.44</td><td align="right">40.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.76</strong></td><td align="right">3.95</td><td align="right">-1.19</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.03</strong></td><td align="right">1.23</td><td align="right">-0.20</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.74</strong></td><td align="right">1.55</td><td align="right">-0.81</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>4.05</strong></td><td align="right">4.60</td><td align="right">-0.55</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right"><strong>2.90</strong></td><td align="right">3.52</td><td align="right">-0.62</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.74</strong></td><td align="right">1.08</td><td align="right">-0.34</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.85</strong></td><td align="right"><strong>2.58</strong></td><td align="right"><strong>-0.74</strong></td><td align="right"><strong>5.71</strong></td><td align="right"><strong>0.00</strong></td></tr>
</tbody>
</table>
