# Distance-ranked Sparse Control

This sanity check compares the real MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.90`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 19/28 instances.**
- **Mean average best deviation: MSA 3.21%, distance-ranked sparse 3.76%, delta -0.54 pp.**
- **Mean success rate: MSA 5.71%, distance-ranked sparse 4.17%, delta +1.55 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.052239.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">4.37</td><td align="right"><strong>4.24</strong></td><td align="right">+0.13</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">93.33</td><td align="right">93.33</td></tr>
<tr><td>crane100_0</td><td align="right">4.52</td><td align="right"><strong>3.79</strong></td><td align="right">+0.73</td><td align="right">0.00</td><td align="right">3.33</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.74</strong></td><td align="right">4.34</td><td align="right">-0.60</td><td align="right">0.00</td><td align="right">3.33</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.70</strong></td><td align="right">3.84</td><td align="right">-1.14</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.44</strong></td><td align="right">4.89</td><td align="right">-0.45</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_1</td><td align="right">3.10</td><td align="right"><strong>2.40</strong></td><td align="right">+0.70</td><td align="right">3.33</td><td align="right">6.67</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.60</strong></td><td align="right">4.43</td><td align="right">-1.83</td><td align="right">13.33</td><td align="right">3.33</td></tr>
<tr><td>dc112</td><td align="right">1.05</td><td align="right"><strong>1.03</strong></td><td align="right">+0.02</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc126</td><td align="right"><strong>2.08</strong></td><td align="right">2.14</td><td align="right">-0.06</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.17</strong></td><td align="right">1.28</td><td align="right">-0.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc176</td><td align="right">2.39</td><td align="right"><strong>2.37</strong></td><td align="right">+0.02</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc188</td><td align="right"><strong>1.51</strong></td><td align="right">1.59</td><td align="right">-0.08</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft53</td><td align="right"><strong>4.90</strong></td><td align="right">6.84</td><td align="right">-1.94</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.24</strong></td><td align="right">2.45</td><td align="right">-0.21</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.62</strong></td><td align="right">6.04</td><td align="right">-2.42</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.72</strong></td><td align="right">6.02</td><td align="right">-0.30</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv120</td><td align="right">4.90</td><td align="right"><strong>4.76</strong></td><td align="right">+0.14</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv130</td><td align="right"><strong>5.12</strong></td><td align="right">5.38</td><td align="right">-0.26</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.03</strong></td><td align="right">6.36</td><td align="right">-1.33</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv150</td><td align="right"><strong>4.52</strong></td><td align="right">5.13</td><td align="right">-0.61</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv160</td><td align="right">4.09</td><td align="right"><strong>3.71</strong></td><td align="right">+0.38</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.94</strong></td><td align="right">5.78</td><td align="right">-0.84</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv55</td><td align="right"><strong>2.32</strong></td><td align="right">2.74</td><td align="right">-0.42</td><td align="right">10.00</td><td align="right">3.33</td></tr>
<tr><td>ftv64</td><td align="right"><strong>1.07</strong></td><td align="right">1.82</td><td align="right">-0.75</td><td align="right">36.67</td><td align="right">0.00</td></tr>
<tr><td>ftv70</td><td align="right">2.92</td><td align="right"><strong>2.50</strong></td><td align="right">+0.42</td><td align="right">3.33</td><td align="right">3.33</td></tr>
<tr><td>ftv90</td><td align="right"><strong>4.04</strong></td><td align="right">8.15</td><td align="right">-4.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.86</strong></td><td align="right">1.16</td><td align="right">-0.30</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>3.21</strong></td><td align="right"><strong>3.76</strong></td><td align="right"><strong>-0.54</strong></td><td align="right"><strong>5.71</strong></td><td align="right"><strong>4.17</strong></td></tr>
</tbody>
</table>
