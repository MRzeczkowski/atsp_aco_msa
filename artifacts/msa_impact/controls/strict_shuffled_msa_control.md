# Strict MSA Shuffled Control

This sanity check compares Strict MSA against deterministic shuffles of the strict MSA mask. Each shuffle preserves the number of strict-MSA boosted directed edges, but assigns them to shuffled directed edges. The control uses the best Strict MSA impact heuristic weight selected separately for each instance.

## Findings

- **Strict MSA had lower average best deviation than the strict shuffled MSA mean in 7/7 instances.**
- **Mean average best deviation: Strict MSA 1.67%, strict shuffled MSA 4.88%, delta -3.21 pp.**
- **Mean success rate: Strict MSA 10.00%, strict shuffled MSA 0.00%, delta +10.00 pp.**
- Strict MSA also beat the best strict shuffled MSA seed in 7/7 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.015625.

## Per-instance comparison

Negative delta means Strict MSA had lower average best deviation than the strict shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>Strict MSA avg best dev. [%]</th><th>strict shuffled MSA mean avg best dev. [%]</th><th>Best strict shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>Strict MSA success [%]</th><th>strict shuffled MSA success [%]</th><th>Seeds</th></tr>
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
