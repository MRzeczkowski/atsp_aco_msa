# MSA Heuristic And Cycle-Cover Overlap

This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.

## Findings

- **34.82% of MSA heuristic edges are also cycle-cover edges.**
- **25.55% of cycle-cover edges are also MSA heuristic edges.**
- **Found-optimal edge partition: both 759, only MSA heuristic 636, only cycle cover 1367, neither 45311.**

<table>
<thead>
<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right">100.00</td><td align="right">43.75</td><td align="right">4</td><td align="right">0</td><td align="right">4</td></tr>
<tr><td>atex3</td><td align="right">100.00</td><td align="right">50.00</td><td align="right">16</td><td align="right">0</td><td align="right">16</td></tr>
<tr><td>atex4</td><td align="right">71.43</td><td align="right">41.67</td><td align="right">20</td><td align="right">8</td><td align="right">28</td></tr>
<tr><td>code198</td><td align="right">1.07</td><td align="right">1.01</td><td align="right">1</td><td align="right">135</td><td align="right">193</td></tr>
<tr><td>crane100_0</td><td align="right">70.91</td><td align="right">39.00</td><td align="right">25</td><td align="right">5</td><td align="right">28</td></tr>
<tr><td>crane100_2</td><td align="right">58.33</td><td align="right">35.00</td><td align="right">29</td><td align="right">9</td><td align="right">29</td></tr>
<tr><td>crane66_0</td><td align="right">58.33</td><td align="right">31.82</td><td align="right">14</td><td align="right">8</td><td align="right">13</td></tr>
<tr><td>crane66_2</td><td align="right">60.53</td><td align="right">34.85</td><td align="right">19</td><td align="right">10</td><td align="right">13</td></tr>
<tr><td>dc112</td><td align="right">2.78</td><td align="right">2.68</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>dc126</td><td align="right">8.33</td><td align="right">7.14</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>dc176</td><td align="right">3.27</td><td align="right">2.84</td><td align="right">0</td><td align="right">0</td><td align="right">0</td></tr>
<tr><td>ft70</td><td align="right">50.94</td><td align="right">38.57</td><td align="right">26</td><td align="right">2</td><td align="right">20</td></tr>
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
<tr><td>ftv55</td><td align="right">77.42</td><td align="right">42.86</td><td align="right">21</td><td align="right">4</td><td align="right">20</td></tr>
<tr><td>ftv70</td><td align="right">76.19</td><td align="right">45.07</td><td align="right">25</td><td align="right">5</td><td align="right">25</td></tr>
<tr><td>rbg358</td><td align="right">22.62</td><td align="right">19.27</td><td align="right">67</td><td align="right">105</td><td align="right">160</td></tr>
<tr><td>rbg403</td><td align="right">16.62</td><td align="right">13.90</td><td align="right">55</td><td align="right">139</td><td align="right">210</td></tr>
<tr><td>rbg443</td><td align="right">20.16</td><td align="right">16.93</td><td align="right">75</td><td align="right">148</td><td align="right">204</td></tr>
<tr><td><strong>Total</strong></td><td align="right">34.82</td><td align="right">25.55</td><td align="right">759</td><td align="right">636</td><td align="right">1367</td></tr>
</tbody>
</table>
