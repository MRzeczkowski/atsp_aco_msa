# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the best MSA-impact heuristic weight selected separately for each instance.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 7/7 instances.**
- **Mean average best deviation: MSA 1.85%, shuffled MSA 4.33%, delta -2.49 pp.**
- **Mean success rate: MSA 5.71%, shuffled MSA 0.00%, delta +5.71 pp.**
- MSA also beat the best shuffled MSA seed in 7/7 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>shuffled MSA mean avg best dev. [%]</th><th>Best shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>0.70</strong></td><td align="right">4.47</td><td align="right">4.29 (seed 102)</td><td align="right">-3.77</td><td align="right">40.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.76</strong></td><td align="right">6.69</td><td align="right">6.10 (seed 101)</td><td align="right">-3.93</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.03</strong></td><td align="right">1.75</td><td align="right">1.46 (seed 101)</td><td align="right">-0.72</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.74</strong></td><td align="right">3.51</td><td align="right">3.16 (seed 101)</td><td align="right">-2.77</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>4.05</strong></td><td align="right">5.08</td><td align="right">4.53 (seed 102)</td><td align="right">-1.03</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ry48p</td><td align="right"><strong>2.90</strong></td><td align="right">6.69</td><td align="right">6.22 (seed 103)</td><td align="right">-3.79</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.74</strong></td><td align="right">2.16</td><td align="right">1.95 (seed 101)</td><td align="right">-1.42</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.85</strong></td><td align="right"><strong>4.33</strong></td><td></td><td align="right"><strong>-2.49</strong></td><td align="right"><strong>5.71</strong></td><td align="right"><strong>0.00</strong></td><td></td></tr>
</tbody>
</table>
