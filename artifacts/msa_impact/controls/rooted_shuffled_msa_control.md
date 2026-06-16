# Rooted MSA Shuffled Control

This sanity check compares Rooted MSA against deterministic per-root shuffles of the rooted MSA masks. Each root preserves the number of boosted directed edges from that root's MSA, but assigns them to shuffled directed edges. The control uses the best Rooted MSA impact heuristic weight selected separately for each instance.

## Findings

- **Rooted MSA had lower average best deviation than the rooted shuffled MSA mean in 7/7 instances.**
- **Mean average best deviation: Rooted MSA 1.65%, rooted shuffled MSA 4.23%, delta -2.58 pp.**
- **Mean success rate: Rooted MSA 7.14%, rooted shuffled MSA 3.33%, delta +3.81 pp.**
- Rooted MSA also beat the best rooted shuffled MSA seed in 7/7 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means Rooted MSA had lower average best deviation than the rooted shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>Rooted MSA avg best dev. [%]</th><th>rooted shuffled MSA mean avg best dev. [%]</th><th>Best rooted shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>Rooted MSA success [%]</th><th>rooted shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>0.76</strong></td><td align="right">3.61</td><td align="right">3.35 (seed 103)</td><td align="right">-2.85</td><td align="right">10.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>2.28</strong></td><td align="right">4.84</td><td align="right">4.07 (seed 101)</td><td align="right">-2.56</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.05</strong></td><td align="right">1.40</td><td align="right">1.35 (seed 101)</td><td align="right">-0.35</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>0.45</strong></td><td align="right">1.23</td><td align="right">0.86 (seed 103)</td><td align="right">-0.78</td><td align="right">40.00</td><td align="right">23.33</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>3.67</strong></td><td align="right">4.23</td><td align="right">4.18 (seed 102)</td><td align="right">-0.56</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ry48p</td><td align="right"><strong>2.67</strong></td><td align="right">11.46</td><td align="right">10.76 (seed 101)</td><td align="right">-8.79</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.70</strong></td><td align="right">2.83</td><td align="right">2.68 (seed 101)</td><td align="right">-2.13</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>1.65</strong></td><td align="right"><strong>4.23</strong></td><td></td><td align="right"><strong>-2.58</strong></td><td align="right"><strong>7.14</strong></td><td align="right"><strong>3.33</strong></td><td></td></tr>
</tbody>
</table>
