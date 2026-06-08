# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.40`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 7/8 instances.**
- **Mean average best deviation: MSA 1.81%, distance-ranked sparse 2.54%, delta -0.73 pp.**
- **Mean success rate: MSA 18.33%, distance-ranked sparse 15.42%, delta +2.92 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">2.78</td><td align="right">-0.66</td><td align="right">3.33</td><td align="right">13.33</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">100.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.51</strong></td><td align="right">5.19</td><td align="right">-1.68</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.10</strong></td><td align="right">1.22</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.56</strong></td><td align="right">0.86</td><td align="right">-0.30</td><td align="right">43.33</td><td align="right">10.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.93</strong></td><td align="right">5.42</td><td align="right">-2.49</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.55</strong></td><td align="right">3.76</td><td align="right">-0.21</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.73</strong></td><td align="right">1.09</td><td align="right">-0.36</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.81</strong></td><td align="right"><strong>2.54</strong></td><td align="right"><strong>-0.73</strong></td><td align="right"><strong>18.33</strong></td><td align="right"><strong>15.42</strong></td></tr>
</tbody>
</table>
