# MSA Impact Summary

Negative deviation delta means the MSA variant had lower average best deviation than the baseline. Strict MSA uses only edges present in every rooted MSA. Rooted MSA gives each ant the MSA rooted at its start vertex. Both variants use the best heuristic weight from the 0.10-1.00 sweep.

## Findings

- **Strict MSA improved average best deviation in 6/7 instances, tied in 1, and was worse in 0.**
- **Mean average best deviation: baseline 2.36%, Strict MSA 1.67%, delta -0.69 pp; success-rate delta +8.57 pp.**
- **Best Strict MSA weights: 0.20 (1), 0.40 (2), 0.50 (1), 0.80 (2), 1.00 (1).**
- **Rooted MSA improved average best deviation in 7/7 instances, tied in 0, and was worse in 0.**
- **Mean average best deviation: baseline 2.36%, Rooted MSA 1.65%, delta -0.71 pp; success-rate delta +5.71 pp.**
- **Best Rooted MSA weights: 0.10 (2), 0.20 (1), 0.30 (1), 0.40 (1), 0.80 (1), 1.00 (1).**

<table>
<thead>
<tr><th>Instance</th><th>Baseline avg best dev. [%]</th><th>Strict MSA avg best dev. [%]</th><th>Strict weight</th><th>Strict delta [pp]</th><th>Rooted MSA avg best dev. [%]</th><th>Rooted weight</th><th>Rooted delta [pp]</th><th>Baseline success [%]</th><th>Strict success [%]</th><th>Rooted success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">2.12</td><td align="right">2.12</td><td align="right">1.00</td><td align="right">+0.00</td><td align="right">0.76</td><td align="right">0.40</td><td align="right">-1.36</td><td align="right">10.00</td><td align="right">20.00</td><td align="right">10.00</td></tr>
<tr><td>ftv64</td><td align="right">1.54</td><td align="right">0.36</td><td align="right">0.50</td><td align="right">-1.18</td><td align="right">0.45</td><td align="right">0.30</td><td align="right">-1.09</td><td align="right">0.00</td><td align="right">50.00</td><td align="right">40.00</td></tr>
<tr><td>ftv90</td><td align="right">3.72</td><td align="right">2.10</td><td align="right">0.80</td><td align="right">-1.62</td><td align="right">3.67</td><td align="right">0.10</td><td align="right">-0.05</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right">3.25</td><td align="right">2.25</td><td align="right">0.40</td><td align="right">-1.00</td><td align="right">2.28</td><td align="right">0.10</td><td align="right">-0.97</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right">0.77</td><td align="right">0.53</td><td align="right">0.80</td><td align="right">-0.24</td><td align="right">0.70</td><td align="right">0.80</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right">3.96</td><td align="right">3.32</td><td align="right">0.40</td><td align="right">-0.64</td><td align="right">2.67</td><td align="right">1.00</td><td align="right">-1.29</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right">1.17</td><td align="right">1.03</td><td align="right">0.20</td><td align="right">-0.14</td><td align="right">1.05</td><td align="right">0.20</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right">2.36</td><td align="right">1.67</td><td></td><td align="right">-0.69</td><td align="right">1.65</td><td></td><td align="right">-0.71</td><td align="right">1.43</td><td align="right">10.00</td><td align="right">7.14</td></tr>
</tbody>
</table>
