# MSA Impact Summary

Negative deviation delta means the MSA heuristic had lower average best deviation than the baseline. MSA rows use the best heuristic weight from the temporary 0.10-1.00 MSA-impact sweep for each instance.

## Findings

- **MSA improved average best deviation in 6/8 instances, tied in 1, and was worse in 1.**
- **Mean average best deviation: baseline 2.16%, MSA 1.66%, delta -0.51 pp.**
- **Mean success-rate delta: +5.83 pp.**
- **Best MSA weights: 0.10 (1), 0.40 (2), 0.50 (2), 0.70 (1), 0.80 (1), 0.90 (1).**

<table>
<thead>
<tr><th>Instance</th><th>Baseline avg best dev. [%]</th><th>MSA avg best dev. [%]</th><th>Best MSA weight</th><th>Delta [pp]</th><th>Baseline success [%]</th><th>MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">2.75</td><td align="right">2.06</td><td align="right">0.70</td><td align="right">-0.69</td><td align="right">3.33</td><td align="right">10.00</td></tr>
<tr><td>ftv64</td><td align="right">1.13</td><td align="right">0.37</td><td align="right">0.50</td><td align="right">-0.76</td><td align="right">13.33</td><td align="right">50.00</td></tr>
<tr><td>ftv90</td><td align="right">3.77</td><td align="right">2.34</td><td align="right">0.80</td><td align="right">-1.43</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_1</td><td align="right">3.29</td><td align="right">3.49</td><td align="right">0.10</td><td align="right">+0.20</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right">0.73</td><td align="right">0.66</td><td align="right">0.90</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ry48p</td><td align="right">4.48</td><td align="right">3.24</td><td align="right">0.40</td><td align="right">-1.24</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right">1.16</td><td align="right">1.09</td><td align="right">0.40</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">0.50</td><td align="right">+0.00</td><td align="right">96.67</td><td align="right">100.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right">2.16</td><td align="right">1.66</td><td></td><td align="right">-0.51</td><td align="right">14.17</td><td align="right">20.00</td></tr>
</tbody>
</table>
