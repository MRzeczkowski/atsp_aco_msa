# Distance-ranked Sparse Control

This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=1.00`.

## Findings

- **MSA had lower average best deviation than the distance-ranked sparse control in 15/26 instances.**
- **Mean average best deviation: MSA 5.57%, distance-ranked sparse 5.83%, delta -0.25 pp.**
- **Mean success rate: MSA 5.38%, distance-ranked sparse 4.54%, delta +0.85 pp.**
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.210040.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">98.00</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>atex4</td><td align="right"><strong>1.94</strong></td><td align="right">2.34</td><td align="right">-0.40</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane100_0</td><td align="right">4.11</td><td align="right"><strong>3.60</strong></td><td align="right">+0.51</td><td align="right">4.00</td><td align="right">4.00</td></tr>
<tr><td>crane100_2</td><td align="right">3.00</td><td align="right"><strong>2.76</strong></td><td align="right">+0.24</td><td align="right">0.00</td><td align="right">2.00</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.51</strong></td><td align="right">4.85</td><td align="right">-0.34</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>3.50</strong></td><td align="right">4.74</td><td align="right">-1.24</td><td align="right">12.00</td><td align="right">0.00</td></tr>
<tr><td>dc112</td><td align="right">1.03</td><td align="right"><strong>1.03</strong></td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc126</td><td align="right">2.21</td><td align="right"><strong>2.04</strong></td><td align="right">+0.17</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>dc176</td><td align="right">2.48</td><td align="right"><strong>2.37</strong></td><td align="right">+0.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.19</strong></td><td align="right">2.74</td><td align="right">-0.55</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.50</strong></td><td align="right">5.14</td><td align="right">-1.64</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.59</strong></td><td align="right">6.08</td><td align="right">-0.49</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.89</strong></td><td align="right">5.06</td><td align="right">-0.17</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.70</strong></td><td align="right">6.56</td><td align="right">-0.86</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.93</strong></td><td align="right">4.69</td><td align="right">-0.76</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv170</td><td align="right"><strong>5.40</strong></td><td align="right">5.84</td><td align="right">-0.44</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv35</td><td align="right"><strong>1.71</strong></td><td align="right">2.26</td><td align="right">-0.55</td><td align="right">8.00</td><td align="right">2.00</td></tr>
<tr><td>ftv38</td><td align="right"><strong>1.84</strong></td><td align="right">3.32</td><td align="right">-1.48</td><td align="right">4.00</td><td align="right">4.00</td></tr>
<tr><td>ftv44</td><td align="right">2.88</td><td align="right"><strong>2.68</strong></td><td align="right">+0.20</td><td align="right">10.00</td><td align="right">4.00</td></tr>
<tr><td>ftv47</td><td align="right"><strong>5.61</strong></td><td align="right">6.20</td><td align="right">-0.59</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv55</td><td align="right"><strong>2.86</strong></td><td align="right">2.92</td><td align="right">-0.06</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>ftv70</td><td align="right">3.58</td><td align="right"><strong>2.38</strong></td><td align="right">+1.20</td><td align="right">2.00</td><td align="right">4.00</td></tr>
<tr><td>rbg358</td><td align="right">19.68</td><td align="right"><strong>19.45</strong></td><td align="right">+0.23</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>rbg403</td><td align="right"><strong>25.24</strong></td><td align="right">25.36</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td>rbg443</td><td align="right">27.26</td><td align="right"><strong>26.79</strong></td><td align="right">+0.47</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>5.57</strong></td><td align="right"><strong>5.83</strong></td><td align="right"><strong>-0.25</strong></td><td align="right"><strong>5.38</strong></td><td align="right"><strong>4.54</strong></td></tr>
</tbody>
</table>
