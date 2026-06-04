# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.90`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 19/28 instances.**
- **Mean average best deviation: MSA 3.22%, distance-ranked sparse 3.73%, delta -0.51 pp.**
- **Mean success rate: MSA 5.64%, distance-ranked sparse 4.14%, delta +1.50 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.052239.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">4.28</td><td align="right"><strong>4.17</strong></td><td align="right">+0.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">94.00</td><td align="right">92.00</td></tr>
<tr><td>crane100_0</td><td align="right">3.98</td><td align="right"><strong>3.86</strong></td><td align="right">+0.12</td><td align="right">4.00</td><td align="right">2.00</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.77</strong></td><td align="right">4.67</td><td align="right">-0.90</td><td align="right">0.00</td><td align="right">2.00</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.63</strong></td><td align="right">3.91</td><td align="right">-1.28</td><td align="right">2.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.48</strong></td><td align="right">4.75</td><td align="right">-0.27</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_1</td><td align="right">2.97</td><td align="right"><strong>2.29</strong></td><td align="right">+0.68</td><td align="right">6.00</td><td align="right">6.00</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.53</strong></td><td align="right">4.50</td><td align="right">-1.97</td><td align="right">10.00</td><td align="right">8.00</td></tr>
<tr><td>dc112</td><td align="right">1.05</td><td align="right"><strong>1.03</strong></td><td align="right">+0.02</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc126</td><td align="right"><strong>2.05</strong></td><td align="right">2.11</td><td align="right">-0.06</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.18</strong></td><td align="right">1.30</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc176</td><td align="right">2.39</td><td align="right"><strong>2.36</strong></td><td align="right">+0.03</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc188</td><td align="right"><strong>1.52</strong></td><td align="right">1.56</td><td align="right">-0.04</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft53</td><td align="right"><strong>5.24</strong></td><td align="right">6.57</td><td align="right">-1.33</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.30</strong></td><td align="right">2.58</td><td align="right">-0.28</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.65</strong></td><td align="right">5.80</td><td align="right">-2.15</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.53</strong></td><td align="right">6.26</td><td align="right">-0.73</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv120</td><td align="right">5.12</td><td align="right"><strong>4.36</strong></td><td align="right">+0.76</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.98</strong></td><td align="right">5.61</td><td align="right">-0.63</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.35</strong></td><td align="right">6.13</td><td align="right">-0.78</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv150</td><td align="right"><strong>4.26</strong></td><td align="right">5.11</td><td align="right">-0.85</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv160</td><td align="right">4.31</td><td align="right"><strong>3.93</strong></td><td align="right">+0.38</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.99</strong></td><td align="right">5.86</td><td align="right">-0.87</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv55</td><td align="right"><strong>2.27</strong></td><td align="right">2.79</td><td align="right">-0.52</td><td align="right">8.00</td><td align="right">4.00</td></tr>
<tr><td>ftv64</td><td align="right"><strong>1.01</strong></td><td align="right">1.72</td><td align="right">-0.71</td><td align="right">32.00</td><td align="right">0.00</td></tr>
<tr><td>ftv70</td><td align="right">3.17</td><td align="right"><strong>2.54</strong></td><td align="right">+0.63</td><td align="right">2.00</td><td align="right">2.00</td></tr>
<tr><td>ftv90</td><td align="right"><strong>4.33</strong></td><td align="right">7.53</td><td align="right">-3.20</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.82</strong></td><td align="right">1.17</td><td align="right">-0.35</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>3.22</strong></td><td align="right"><strong>3.73</strong></td><td align="right"><strong>-0.51</strong></td><td align="right"><strong>5.64</strong></td><td align="right"><strong>4.14</strong></td></tr>
</tbody>
</table>
