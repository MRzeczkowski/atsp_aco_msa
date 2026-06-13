# MSA Impact Summary

Negative deviation delta means the MSA heuristic had lower average best deviation than the baseline. In this temporary pipeline the MSA heuristic uses the thinnest cached root MSA for each instance and the best heuristic weight from the 0.10-1.00 sweep.

## Findings

- **MSA heuristic improved average best deviation in 6/7 instances, tied in 0, and was worse in 1.**
- **Mean average best deviation: baseline 2.36%, MSA heuristic 1.85%, delta -0.52 pp.**
- **Mean success-rate delta: +4.29 pp.**
- **Best MSA heuristic weights: 0.20 (2), 0.30 (3), 0.50 (1), 0.80 (1).**

<table>
<thead>
<tr><th>Instance</th><th>Baseline avg best dev. [%]</th><th>MSA heuristic avg best dev. [%]</th><th>Best MSA weight</th><th>Delta [pp]</th><th>Baseline success [%]</th><th>MSA heuristic success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">2.12</td><td align="right">0.70</td><td align="right">0.50</td><td align="right">-1.42</td><td align="right">10.00</td><td align="right">40.00</td></tr>
<tr><td>ftv64</td><td align="right">1.54</td><td align="right">0.74</td><td align="right">0.80</td><td align="right">-0.80</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv90</td><td align="right">3.72</td><td align="right">4.05</td><td align="right">0.30</td><td align="right">+0.33</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right">3.25</td><td align="right">2.76</td><td align="right">0.20</td><td align="right">-0.49</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right">0.77</td><td align="right">0.74</td><td align="right">0.30</td><td align="right">-0.03</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right">3.96</td><td align="right">2.90</td><td align="right">0.20</td><td align="right">-1.06</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right">1.17</td><td align="right">1.03</td><td align="right">0.30</td><td align="right">-0.14</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right">2.36</td><td align="right">1.85</td><td></td><td align="right">-0.52</td><td align="right">1.43</td><td align="right">5.71</td></tr>
</tbody>
</table>
