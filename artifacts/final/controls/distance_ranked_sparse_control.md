# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.40`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 16/27 instances.**
- **Mean average best deviation: MSA 4.46%, distance-ranked sparse 4.64%, delta -0.17 pp.**
- **Mean success rate: MSA 9.78%, distance-ranked sparse 8.96%, delta +0.81 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.151590.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">100.00</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>atex4</td><td align="right"><strong>2.00</strong></td><td align="right">2.07</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">100.00</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.91</strong></td><td align="right">4.16</td><td align="right">-0.25</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.79</strong></td><td align="right">2.88</td><td align="right">-0.09</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_0</td><td align="right">2.88</td><td align="right"><strong>2.83</strong></td><td align="right">+0.05</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.05</strong></td><td align="right">2.50</td><td align="right">-0.45</td><td align="right">36.00</td><td align="right">22.00</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.01</strong></td><td align="right">1.03</td><td align="right">-0.02</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.84</strong></td><td align="right">1.85</td><td align="right">-0.01</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.12</strong></td><td align="right">2.25</td><td align="right">-0.13</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.50</strong></td><td align="right">1.96</td><td align="right">-0.46</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv100</td><td align="right"><strong>2.70</strong></td><td align="right">5.22</td><td align="right">-2.52</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.19</strong></td><td align="right">6.16</td><td align="right">-0.97</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.40</strong></td><td align="right">5.19</td><td align="right">-0.79</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv140</td><td align="right">5.31</td><td align="right"><strong>5.29</strong></td><td align="right">+0.02</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.36</strong></td><td align="right">3.53</td><td align="right">-0.17</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv170</td><td align="right">5.11</td><td align="right"><strong>4.59</strong></td><td align="right">+0.52</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv35</td><td align="right"><strong>2.06</strong></td><td align="right">2.30</td><td align="right">-0.24</td><td align="right">4.00</td><td align="right">2.00</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.14</strong></td><td align="right">3.11</td><td align="right">-0.97</td><td align="right">2.00</td><td align="right">8.00</td></tr>
<tr><td>ftv44</td><td align="right">2.89</td><td align="right"><strong>2.41</strong></td><td align="right">+0.48</td><td align="right">4.00</td><td align="right">4.00</td></tr>
<tr><td>ftv47</td><td align="right"><strong>5.51</strong></td><td align="right">5.84</td><td align="right">-0.33</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.91</strong></td><td align="right">2.26</td><td align="right">-0.35</td><td align="right">8.00</td><td align="right">2.00</td></tr>
<tr><td>ftv70</td><td align="right">2.17</td><td align="right"><strong>1.39</strong></td><td align="right">+0.78</td><td align="right">10.00</td><td align="right">4.00</td></tr>
<tr><td>rbg358</td><td align="right">14.98</td><td align="right"><strong>14.36</strong></td><td align="right">+0.62</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>rbg403</td><td align="right">21.28</td><td align="right"><strong>21.20</strong></td><td align="right">+0.08</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>rbg443</td><td align="right">21.15</td><td align="right"><strong>20.53</strong></td><td align="right">+0.62</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.46</strong></td><td align="right"><strong>4.64</strong></td><td align="right"><strong>-0.17</strong></td><td align="right"><strong>9.78</strong></td><td align="right"><strong>8.96</strong></td></tr>
</tbody>
</table>
