# MSA Heuristic And Cycle-Cover Overlap

This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.

## Findings

- **35.48% of MSA heuristic edges are also cycle-cover edges.**
- **22.54% of cycle-cover edges are also MSA heuristic edges.**
- **Found-optimal edge partition: both 307, only MSA heuristic 143, only cycle cover 606, neither 27843.**

<table>
<thead>
<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">100.00</td><td align="right">44.44</td><td align="right">32</td><td align="right">0</td><td align="right">40</td></tr>
<tr><td>code198</td><td align="right">1.07</td><td align="right">1.01</td><td align="right">1</td><td align="right">70</td><td align="right">193</td></tr>
<tr><td>crane100_1</td><td align="right">56.14</td><td align="right">32.00</td><td align="right">27</td><td align="right">6</td><td align="right">33</td></tr>
<tr><td>crane66_1</td><td align="right">60.61</td><td align="right">30.30</td><td align="right">16</td><td align="right">2</td><td align="right">20</td></tr>
<tr><td>dc134</td><td align="right">8.22</td><td align="right">4.48</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>dc188</td><td align="right">2.19</td><td align="right">2.13</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>ft53</td><td align="right">46.43</td><td align="right">24.53</td><td align="right">9</td><td align="right">6</td><td align="right">25</td></tr>
<tr><td>ftv120</td><td align="right">71.23</td><td align="right">42.98</td><td align="right">45</td><td align="right">5</td><td align="right">48</td></tr>
<tr><td>ftv150</td><td align="right">71.74</td><td align="right">43.71</td><td align="right">61</td><td align="right">9</td><td align="right">66</td></tr>
<tr><td>ftv33</td><td align="right">70.00</td><td align="right">41.18</td><td align="right">11</td><td align="right">3</td><td align="right">7</td></tr>
<tr><td>ftv64</td><td align="right">71.05</td><td align="right">41.54</td><td align="right">25</td><td align="right">4</td><td align="right">23</td></tr>
<tr><td>ftv90</td><td align="right">75.00</td><td align="right">42.86</td><td align="right">35</td><td align="right">3</td><td align="right">40</td></tr>
<tr><td>p43</td><td align="right">63.64</td><td align="right">32.56</td><td align="right">14</td><td align="right">8</td><td align="right">24</td></tr>
<tr><td>rbg323</td><td align="right">31.29</td><td align="right">14.24</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>ry48p</td><td align="right">64.71</td><td align="right">22.92</td><td align="right">6</td><td align="right">5</td><td align="right">14</td></tr>
<tr><td>td100_1</td><td align="right">30.49</td><td align="right">24.75</td><td align="right">25</td><td align="right">22</td><td align="right">73</td></tr>
<tr><td><strong>Total</strong></td><td align="right">35.48</td><td align="right">22.54</td><td align="right">307</td><td align="right">143</td><td align="right">606</td></tr>
</tbody>
</table>
