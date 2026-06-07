# MSA Impact Structure

This report describes the active MSA heuristic modifier matrix used by the MSA impact pipeline. Negative deviation delta means the MSA heuristic had lower average best deviation than the baseline.

## Findings

- **Average boosted-edge ratio against n-1: 0.61.**
- **Average missing outgoing vertices: 75.38; average missing incoming vertices: 33.88.**
- **Average missing endpoints: improved cases 94.80, tied/worse cases 133.33.**

<table>
<thead>
<tr><th>Instance</th><th>n</th><th>Boosted edges</th><th>Boosted/(n-1)</th><th>Missing outgoing</th><th>Missing incoming</th><th>Dev. delta [pp]</th><th>Success delta [pp]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">72</td><td align="right">32</td><td align="right">0.45</td><td align="right">40</td><td align="right">40</td><td align="right">-0.63</td><td align="right">+0.00</td></tr>
<tr><td>ftv64</td><td align="right">65</td><td align="right">38</td><td align="right">0.59</td><td align="right">34</td><td align="right">27</td><td align="right">-0.14</td><td align="right">+26.67</td></tr>
<tr><td>ftv90</td><td align="right">91</td><td align="right">52</td><td align="right">0.58</td><td align="right">45</td><td align="right">39</td><td align="right">-1.29</td><td align="right">+0.00</td></tr>
<tr><td>crane100_1</td><td align="right">100</td><td align="right">57</td><td align="right">0.58</td><td align="right">59</td><td align="right">43</td><td align="right">+0.48</td><td align="right">+0.00</td></tr>
<tr><td>td100_1</td><td align="right">101</td><td align="right">82</td><td align="right">0.82</td><td align="right">72</td><td align="right">19</td><td align="right">+0.01</td><td align="right">+0.00</td></tr>
<tr><td>ry48p</td><td align="right">48</td><td align="right">17</td><td align="right">0.36</td><td align="right">32</td><td align="right">31</td><td align="right">-1.24</td><td align="right">+0.00</td></tr>
<tr><td>dc134</td><td align="right">134</td><td align="right">73</td><td align="right">0.55</td><td align="right">125</td><td align="right">61</td><td align="right">-0.07</td><td align="right">+0.00</td></tr>
<tr><td>code198</td><td align="right">198</td><td align="right">187</td><td align="right">0.95</td><td align="right">196</td><td align="right">11</td><td align="right">+0.00</td><td align="right">+3.33</td></tr>
</tbody>
</table>
