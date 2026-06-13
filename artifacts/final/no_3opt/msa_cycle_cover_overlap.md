# MSA Heuristic And Cycle-Cover Overlap

This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.

## Findings

- **42.41% of MSA heuristic edges are also cycle-cover edges.**
- **28.81% of cycle-cover edges are also MSA heuristic edges.**
- **Found-optimal edge partition: both 151, only MSA heuristic 41, only cycle cover 222, neither 579.**

<table>
<thead>
<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">100.00</td><td align="right">44.44</td><td align="right">32</td><td align="right">0</td><td align="right">40</td></tr>
<tr><td>crane100_1</td><td align="right">56.25</td><td align="right">36.00</td><td align="right">28</td><td align="right">7</td><td align="right">32</td></tr>
<tr><td>dc134</td><td align="right">4.65</td><td align="right">4.48</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>ftv64</td><td align="right">71.05</td><td align="right">41.54</td><td align="right">25</td><td align="right">4</td><td align="right">23</td></tr>
<tr><td>ftv90</td><td align="right">75.00</td><td align="right">42.86</td><td align="right">35</td><td align="right">3</td><td align="right">40</td></tr>
<tr><td>ry48p</td><td align="right">64.71</td><td align="right">22.92</td><td align="right">6</td><td align="right">5</td><td align="right">14</td></tr>
<tr><td>td100_1</td><td align="right">30.12</td><td align="right">24.75</td><td align="right">25</td><td align="right">22</td><td align="right">73</td></tr>
<tr><td><strong>Total</strong></td><td align="right">42.41</td><td align="right">28.81</td><td align="right">151</td><td align="right">41</td><td align="right">222</td></tr>
</tbody>
</table>
