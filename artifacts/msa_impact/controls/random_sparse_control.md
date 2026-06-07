# Random Sparse Control

This sanity check compares the final MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison reads MSA from the final results and averages the available final-control random seeds for each instance.

## Findings

- **MSA had lower average best deviation than the random-sparse mean in 8/8 instances.**
- **Mean average best deviation: MSA 1.80%, random sparse 3.60%, delta -1.80 pp.**
- **Mean success rate: MSA 17.92%, random sparse 10.97%, delta +6.94 pp.**
- MSA also beat the best random seed in 6/8 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.007812.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">2.66</td><td align="right">1.88 (seed 3)</td><td align="right">-0.54</td><td align="right">3.33</td><td align="right">1.11</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 1)</td><td align="right">-0.00</td><td align="right">100.00</td><td align="right">81.11</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.77</strong></td><td align="right">8.25</td><td align="right">7.15 (seed 3)</td><td align="right">-4.48</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.09</strong></td><td align="right">1.90</td><td align="right">1.79 (seed 2)</td><td align="right">-0.81</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.99</strong></td><td align="right">1.75</td><td align="right">1.52 (seed 1)</td><td align="right">-0.76</td><td align="right">40.00</td><td align="right">5.56</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.48</strong></td><td align="right">5.90</td><td align="right">4.92 (seed 1)</td><td align="right">-3.42</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.24</strong></td><td align="right">6.31</td><td align="right">5.88 (seed 2)</td><td align="right">-3.07</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.74</strong></td><td align="right">2.07</td><td align="right">1.89 (seed 1)</td><td align="right">-1.33</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.80</strong></td><td align="right"><strong>3.60</strong></td><td></td><td align="right"><strong>-1.80</strong></td><td align="right"><strong>17.92</strong></td><td align="right"><strong>10.97</strong></td><td></td></tr>
</tbody>
</table>
