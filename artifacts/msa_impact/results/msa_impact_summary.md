# MSA Impact Summary

Negative deviation delta means the MSA heuristic had lower average best deviation than the baseline. In this temporary pipeline each ant uses the MSA rooted at its start vertex and the best heuristic weight from the 0.10-1.00 sweep.

## Findings

- **MSA heuristic improved average best deviation in 7/7 instances, tied in 0, and was worse in 0.**
- **Mean average best deviation: baseline 2.36%, MSA heuristic 1.65%, delta -0.71 pp.**
- **Mean success-rate delta: +5.71 pp.**
- **Best MSA heuristic weights: 0.10 (2), 0.20 (1), 0.30 (1), 0.40 (1), 0.80 (1), 1.00 (1).**

<table>
<thead>
<tr><th>Instance</th><th>Baseline avg best dev. [%]</th><th>MSA heuristic avg best dev. [%]</th><th>Best MSA weight</th><th>Delta [pp]</th><th>Baseline success [%]</th><th>MSA heuristic success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">2.12</td><td align="right">0.76</td><td align="right">0.40</td><td align="right">-1.36</td><td align="right">10.00</td><td align="right">10.00</td></tr>
<tr><td>ftv64</td><td align="right">1.54</td><td align="right">0.45</td><td align="right">0.30</td><td align="right">-1.09</td><td align="right">0.00</td><td align="right">40.00</td></tr>
<tr><td>ftv90</td><td align="right">3.72</td><td align="right">3.67</td><td align="right">0.10</td><td align="right">-0.05</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right">3.25</td><td align="right">2.28</td><td align="right">0.10</td><td align="right">-0.97</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right">0.77</td><td align="right">0.70</td><td align="right">0.80</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right">3.96</td><td align="right">2.67</td><td align="right">1.00</td><td align="right">-1.29</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right">1.17</td><td align="right">1.05</td><td align="right">0.20</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right">2.36</td><td align="right">1.65</td><td></td><td align="right">-0.71</td><td align="right">1.43</td><td align="right">7.14</td></tr>
</tbody>
</table>
