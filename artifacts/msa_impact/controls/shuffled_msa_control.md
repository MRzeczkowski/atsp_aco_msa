# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the same `heuristicWeight=0.40`.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 8/8 instances.**
- **Mean average best deviation: MSA 1.80%, shuffled MSA 3.66%, delta -1.86 pp.**
- **Mean success rate: MSA 17.92%, shuffled MSA 11.25%, delta +6.67 pp.**
- MSA also beat the best shuffled MSA seed in 7/8 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.007812.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>shuffled MSA mean avg best dev. [%]</th><th>Best shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>2.12</strong></td><td align="right">3.68</td><td align="right">3.05 (seed 102)</td><td align="right">-1.56</td><td align="right">3.33</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 103)</td><td align="right">-0.00</td><td align="right">100.00</td><td align="right">86.67</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.77</strong></td><td align="right">7.97</td><td align="right">7.19 (seed 103)</td><td align="right">-4.20</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.09</strong></td><td align="right">1.65</td><td align="right">1.51 (seed 101)</td><td align="right">-0.56</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.99</strong></td><td align="right">1.54</td><td align="right">1.36 (seed 103)</td><td align="right">-0.55</td><td align="right">40.00</td><td align="right">3.33</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.48</strong></td><td align="right">5.73</td><td align="right">4.83 (seed 101)</td><td align="right">-3.25</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.24</strong></td><td align="right">6.71</td><td align="right">5.80 (seed 101)</td><td align="right">-3.47</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.74</strong></td><td align="right">1.99</td><td align="right">1.93 (seed 101)</td><td align="right">-1.25</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.80</strong></td><td align="right"><strong>3.66</strong></td><td></td><td align="right"><strong>-1.86</strong></td><td align="right"><strong>17.92</strong></td><td align="right"><strong>11.25</strong></td><td></td></tr>
</tbody>
</table>
