# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.70`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 15/26 instances.**
- **Mean average best deviation: MSA 4.60%, distance-ranked sparse 4.92%, delta -0.32 pp.**
- **Mean success rate: MSA 7.69%, distance-ranked sparse 5.69%, delta +2.00 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.210040.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>0.00</strong></td><td align="right">0.01</td><td align="right">-0.01</td><td align="right">98.00</td><td align="right">96.00</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>atex4</td><td align="right"><strong>1.97</strong></td><td align="right">2.19</td><td align="right">-0.22</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.62</strong></td><td align="right">4.25</td><td align="right">-0.63</td><td align="right">2.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>3.09</strong></td><td align="right">3.19</td><td align="right">-0.10</td><td align="right">2.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.28</strong></td><td align="right">2.90</td><td align="right">-0.62</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.05</strong></td><td align="right">2.86</td><td align="right">-0.81</td><td align="right">42.00</td><td align="right">22.00</td></tr>
<tr><td>dc112</td><td align="right">1.02</td><td align="right"><strong>1.01</strong></td><td align="right">+0.01</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc126</td><td align="right">1.94</td><td align="right"><strong>1.83</strong></td><td align="right">+0.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc176</td><td align="right">2.28</td><td align="right"><strong>2.25</strong></td><td align="right">+0.03</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.56</strong></td><td align="right">1.65</td><td align="right">-0.09</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.40</strong></td><td align="right">6.60</td><td align="right">-3.20</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv110</td><td align="right"><strong>4.42</strong></td><td align="right">5.54</td><td align="right">-1.12</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.44</strong></td><td align="right">5.00</td><td align="right">-0.56</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv140</td><td align="right"><strong>4.85</strong></td><td align="right">5.35</td><td align="right">-0.50</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.15</strong></td><td align="right">3.71</td><td align="right">-0.56</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.83</strong></td><td align="right">5.67</td><td align="right">-0.84</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv35</td><td align="right">2.20</td><td align="right"><strong>2.20</strong></td><td align="right">+0.00</td><td align="right">4.00</td><td align="right">4.00</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.36</strong></td><td align="right">3.04</td><td align="right">-0.68</td><td align="right">2.00</td><td align="right">6.00</td></tr>
<tr><td>ftv44</td><td align="right">2.88</td><td align="right"><strong>2.88</strong></td><td align="right">+0.00</td><td align="right">12.00</td><td align="right">4.00</td></tr>
<tr><td>ftv47</td><td align="right">5.58</td><td align="right"><strong>5.45</strong></td><td align="right">+0.13</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.50</strong></td><td align="right">2.28</td><td align="right">-0.78</td><td align="right">20.00</td><td align="right">6.00</td></tr>
<tr><td>ftv70</td><td align="right">2.40</td><td align="right"><strong>1.46</strong></td><td align="right">+0.94</td><td align="right">18.00</td><td align="right">10.00</td></tr>
<tr><td>rbg358</td><td align="right">14.84</td><td align="right"><strong>14.47</strong></td><td align="right">+0.37</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>rbg403</td><td align="right">21.51</td><td align="right"><strong>21.36</strong></td><td align="right">+0.15</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>rbg443</td><td align="right">21.07</td><td align="right"><strong>20.53</strong></td><td align="right">+0.54</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.60</strong></td><td align="right"><strong>4.92</strong></td><td align="right"><strong>-0.32</strong></td><td align="right"><strong>7.69</strong></td><td align="right"><strong>5.69</strong></td></tr>
</tbody>
</table>
