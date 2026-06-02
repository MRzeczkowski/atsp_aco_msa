# MSA Heuristic And Cycle-Cover Overlap

This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.

## Findings

- **73.06% of MSA heuristic edges are also cycle-cover edges.**
- **39.17% of cycle-cover edges are also MSA heuristic edges.**
- **Found-optimal edge partition: both 127, only MSA heuristic 18, only cycle cover 136, neither 410.**

<table>
<thead>
<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">100.00</td><td align="right">44.44</td><td align="right">32</td><td align="right">0</td><td align="right">40</td></tr>
<tr><td>crane66_1</td><td align="right">60.61</td><td align="right">30.30</td><td align="right">16</td><td align="right">2</td><td align="right">20</td></tr>
<tr><td>crane66_2</td><td align="right">60.53</td><td align="right">34.85</td><td align="right">19</td><td align="right">9</td><td align="right">13</td></tr>
<tr><td>ftv64</td><td align="right">71.05</td><td align="right">41.54</td><td align="right">25</td><td align="right">4</td><td align="right">23</td></tr>
<tr><td>ftv90</td><td align="right">75.00</td><td align="right">42.86</td><td align="right">35</td><td align="right">3</td><td align="right">40</td></tr>
<tr><td><strong>Total</strong></td><td align="right">73.06</td><td align="right">39.17</td><td align="right">127</td><td align="right">18</td><td align="right">136</td></tr>
</tbody>
</table>
