# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the best MSA-impact heuristic weight selected separately for each instance.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 7/7 instances.**
- **Mean average best deviation: MSA 1.67%, shuffled MSA 4.88%, delta -3.21 pp.**
- **Mean success rate: MSA 10.00%, shuffled MSA 0.00%, delta +10.00 pp.**
- MSA also beat the best shuffled MSA seed in 7/7 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>shuffled MSA mean avg best dev. [%]</th><th>Best shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">4.74</td><td align="right">4.56 (seed 102)</td><td align="right">-2.62</td><td align="right">20.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.25</strong></td><td align="right">8.57</td><td align="right">7.89 (seed 103)</td><td align="right">-6.32</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.03</strong></td><td align="right">1.49</td><td align="right">1.41 (seed 101)</td><td align="right">-0.46</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.36</strong></td><td align="right">2.21</td><td align="right">1.52 (seed 102)</td><td align="right">-1.85</td><td align="right">50.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.10</strong></td><td align="right">8.23</td><td align="right">7.62 (seed 102)</td><td align="right">-6.13</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.32</strong></td><td align="right">6.29</td><td align="right">5.93 (seed 101)</td><td align="right">-2.97</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.53</strong></td><td align="right">2.63</td><td align="right">2.38 (seed 103)</td><td align="right">-2.10</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.67</strong></td><td align="right"><strong>4.88</strong></td><td></td><td align="right"><strong>-3.21</strong></td><td align="right"><strong>10.00</strong></td><td align="right"><strong>0.00</strong></td><td></td></tr>
</tbody>
</table>
