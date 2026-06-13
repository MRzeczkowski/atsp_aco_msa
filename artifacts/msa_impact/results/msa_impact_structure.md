# MSA Impact Structure

This report describes the default strict composite MSA modifier matrix used by the MSA impact pipeline. Negative deviation delta means the MSA heuristic had lower average best deviation than the baseline.

## Findings

- **Average boosted-edge ratio against n-1: 0.63.**
- **Average missing outgoing vertices: 57.71; average missing incoming vertices: 28.00.**
- **Average missing endpoints: improved cases 86.67, tied/worse cases 80.00.**

<table>
<thead>
<tr><th>Instance</th><th>n</th><th>Boosted edges</th><th>Boosted/(n-1)</th><th>Missing outgoing</th><th>Missing incoming</th><th>Dev. delta [pp]</th><th>Success delta [pp]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">72</td><td align="right">32</td><td align="right">0.45</td><td align="right">40</td><td align="right">40</td><td align="right">+0.00</td><td align="right">+10.00</td></tr>
<tr><td>ftv64</td><td align="right">65</td><td align="right">38</td><td align="right">0.59</td><td align="right">34</td><td align="right">27</td><td align="right">-1.18</td><td align="right">+60.00</td></tr>
<tr><td>ftv90</td><td align="right">91</td><td align="right">52</td><td align="right">0.58</td><td align="right">45</td><td align="right">39</td><td align="right">-1.64</td><td align="right">+0.00</td></tr>
<tr><td>crane100_1</td><td align="right">100</td><td align="right">64</td><td align="right">0.65</td><td align="right">56</td><td align="right">36</td><td align="right">-0.64</td><td align="right">+0.00</td></tr>
<tr><td>td100_1</td><td align="right">101</td><td align="right">83</td><td align="right">0.83</td><td align="right">72</td><td align="right">18</td><td align="right">-0.12</td><td align="right">+0.00</td></tr>
<tr><td>ry48p</td><td align="right">48</td><td align="right">17</td><td align="right">0.36</td><td align="right">32</td><td align="right">31</td><td align="right">-0.90</td><td align="right">+0.00</td></tr>
<tr><td>dc134</td><td align="right">134</td><td align="right">129</td><td align="right">0.97</td><td align="right">125</td><td align="right">5</td><td align="right">-0.11</td><td align="right">+0.00</td></tr>
</tbody>
</table>
