# MSA Heuristic And Cycle-Cover Overlap

This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.

## Findings

- **37.95% of MSA heuristic edges are also cycle-cover edges.**
- **26.93% of cycle-cover edges are also MSA heuristic edges.**
- **Found-optimal edge partition: both 753, only MSA heuristic 499, only cycle cover 1179, neither 18101.**

<table>
<thead>
<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right">100.00</td><td align="right">43.75</td><td align="right">4</td><td align="right">0</td><td align="right">4</td></tr>
<tr><td>atex3</td><td align="right">100.00</td><td align="right">37.50</td><td align="right">12</td><td align="right">0</td><td align="right">20</td></tr>
<tr><td>atex4</td><td align="right">74.07</td><td align="right">41.67</td><td align="right">20</td><td align="right">7</td><td align="right">28</td></tr>
<tr><td>crane100_0</td><td align="right">69.64</td><td align="right">39.00</td><td align="right">25</td><td align="right">5</td><td align="right">28</td></tr>
<tr><td>crane100_2</td><td align="right">58.33</td><td align="right">35.00</td><td align="right">29</td><td align="right">9</td><td align="right">29</td></tr>
<tr><td>crane66_0</td><td align="right">58.33</td><td align="right">31.82</td><td align="right">14</td><td align="right">8</td><td align="right">13</td></tr>
<tr><td>crane66_2</td><td align="right">60.53</td><td align="right">34.85</td><td align="right">19</td><td align="right">9</td><td align="right">13</td></tr>
<tr><td>dc112</td><td align="right">2.80</td><td align="right">2.68</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>dc126</td><td align="right">8.33</td><td align="right">7.14</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>dc176</td><td align="right">4.10</td><td align="right">2.84</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>ft70</td><td align="right">51.92</td><td align="right">38.57</td><td align="right">26</td><td align="right">1</td><td align="right">20</td></tr>
<tr><td>ftv100</td><td align="right">74.14</td><td align="right">42.57</td><td align="right">39</td><td align="right">4</td><td align="right">48</td></tr>
<tr><td>ftv110</td><td align="right">70.31</td><td align="right">40.54</td><td align="right">38</td><td align="right">6</td><td align="right">50</td></tr>
<tr><td>ftv130</td><td align="right">73.75</td><td align="right">45.04</td><td align="right">54</td><td align="right">6</td><td align="right">52</td></tr>
<tr><td>ftv140</td><td align="right">74.12</td><td align="right">44.68</td><td align="right">57</td><td align="right">7</td><td align="right">57</td></tr>
<tr><td>ftv160</td><td align="right">74.47</td><td align="right">43.48</td><td align="right">65</td><td align="right">7</td><td align="right">78</td></tr>
<tr><td>ftv170</td><td align="right">73.27</td><td align="right">43.27</td><td align="right">67</td><td align="right">10</td><td align="right">74</td></tr>
<tr><td>ftv35</td><td align="right">61.90</td><td align="right">36.11</td><td align="right">9</td><td align="right">4</td><td align="right">7</td></tr>
<tr><td>ftv38</td><td align="right">56.52</td><td align="right">33.33</td><td align="right">9</td><td align="right">5</td><td align="right">9</td></tr>
<tr><td>ftv44</td><td align="right">65.38</td><td align="right">37.78</td><td align="right">13</td><td align="right">5</td><td align="right">10</td></tr>
<tr><td>ftv47</td><td align="right">60.00</td><td align="right">31.25</td><td align="right">11</td><td align="right">4</td><td align="right">19</td></tr>
<tr><td>ftv55</td><td align="right">76.67</td><td align="right">41.07</td><td align="right">20</td><td align="right">4</td><td align="right">21</td></tr>
<tr><td>ftv70</td><td align="right">76.19</td><td align="right">45.07</td><td align="right">25</td><td align="right">5</td><td align="right">25</td></tr>
<tr><td>rbg358</td><td align="right">22.62</td><td align="right">19.27</td><td align="right">66</td><td align="right">105</td><td align="right">161</td></tr>
<tr><td>rbg403</td><td align="right">16.62</td><td align="right">13.90</td><td align="right">55</td><td align="right">139</td><td align="right">210</td></tr>
<tr><td>rbg443</td><td align="right">20.32</td><td align="right">17.16</td><td align="right">76</td><td align="right">149</td><td align="right">203</td></tr>
<tr><td><strong>Total</strong></td><td align="right">37.95</td><td align="right">26.93</td><td align="right">753</td><td align="right">499</td><td align="right">1179</td></tr>
</tbody>
</table>
