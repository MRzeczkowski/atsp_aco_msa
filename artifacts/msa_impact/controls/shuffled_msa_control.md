# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the best MSA-impact heuristic weight selected separately for each instance.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 7/7 instances.**
- **Mean average best deviation: MSA 1.71%, shuffled MSA 4.96%, delta -3.26 pp.**
- **Mean success rate: MSA 11.43%, shuffled MSA 0.48%, delta +10.95 pp.**
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
<tr><td>crane100_1</td><td align="right"><strong>2.61</strong></td><td align="right">4.41</td><td align="right">3.83 (seed 103)</td><td align="right">-1.80</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.06</strong></td><td align="right">1.47</td><td align="right">1.41 (seed 101)</td><td align="right">-0.41</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.36</strong></td><td align="right">2.85</td><td align="right">2.33 (seed 101)</td><td align="right">-2.49</td><td align="right">60.00</td><td align="right">3.33</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>2.08</strong></td><td align="right">9.05</td><td align="right">8.50 (seed 101)</td><td align="right">-6.97</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ry48p</td><td align="right"><strong>3.06</strong></td><td align="right">9.66</td><td align="right">9.14 (seed 102)</td><td align="right">-6.60</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.65</strong></td><td align="right">2.57</td><td align="right">2.54 (seed 102)</td><td align="right">-1.92</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.71</strong></td><td align="right"><strong>4.96</strong></td><td></td><td align="right"><strong>-3.26</strong></td><td align="right"><strong>11.43</strong></td><td align="right"><strong>0.48</strong></td><td></td></tr>
</tbody>
</table>
