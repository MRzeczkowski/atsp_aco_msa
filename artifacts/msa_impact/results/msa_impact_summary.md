# MSA Impact Summary

Negative deviation delta means the MSA heuristic had lower average best deviation than the baseline. In this temporary pipeline the MSA heuristic uses the default strict composite MSA matrix and the best heuristic weight from the 0.10-1.00 sweep.

## Findings

- **MSA heuristic improved average best deviation in 6/7 instances, tied in 1, and was worse in 0.**
- **Mean average best deviation: baseline 2.36%, MSA heuristic 1.71%, delta -0.66 pp.**
- **Mean success-rate delta: +10.00 pp.**
- **Best MSA heuristic weights: 0.10 (1), 0.20 (1), 0.70 (1), 0.90 (2), 1.00 (2).**

<table>
<thead>
<tr><th>Instance</th><th>Baseline avg best dev. [%]</th><th>MSA heuristic avg best dev. [%]</th><th>Best MSA weight</th><th>Delta [pp]</th><th>Baseline success [%]</th><th>MSA heuristic success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">2.12</td><td align="right">2.12</td><td align="right">1.00</td><td align="right">+0.00</td><td align="right">10.00</td><td align="right">20.00</td></tr>
<tr><td>ftv64</td><td align="right">1.54</td><td align="right">0.36</td><td align="right">0.70</td><td align="right">-1.18</td><td align="right">0.00</td><td align="right">60.00</td></tr>
<tr><td>ftv90</td><td align="right">3.72</td><td align="right">2.08</td><td align="right">1.00</td><td align="right">-1.64</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right">3.25</td><td align="right">2.61</td><td align="right">0.10</td><td align="right">-0.64</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right">0.77</td><td align="right">0.65</td><td align="right">0.90</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right">3.96</td><td align="right">3.06</td><td align="right">0.90</td><td align="right">-0.90</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right">1.17</td><td align="right">1.06</td><td align="right">0.20</td><td align="right">-0.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right">2.36</td><td align="right">1.71</td><td></td><td align="right">-0.66</td><td align="right">1.43</td><td align="right">11.43</td></tr>
</tbody>
</table>
