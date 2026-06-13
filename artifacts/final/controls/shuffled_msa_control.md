# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the same `heuristicWeight=0.40`.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 25/27 instances.**
- **Mean average best deviation: MSA 4.46%, shuffled MSA 6.69%, delta -2.22 pp.**
- **Mean success rate: MSA 9.78%, shuffled MSA 7.11%, delta +2.67 pp.**
- MSA also beat the best shuffled MSA seed in 25/27 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.000000.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>shuffled MSA mean avg best dev. [%]</th><th>Best shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>0.00</strong></td><td align="right">0.01</td><td align="right">0.01 (seed 103)</td><td align="right">-0.01</td><td align="right">100.00</td><td align="right">94.00</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 102)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>2.00</strong></td><td align="right">3.51</td><td align="right">2.71 (seed 102)</td><td align="right">-1.51</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">0.00 (seed 103)</td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">88.67</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.91</strong></td><td align="right">6.98</td><td align="right">6.71 (seed 103)</td><td align="right">-3.07</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.79</strong></td><td align="right">6.39</td><td align="right">6.06 (seed 101)</td><td align="right">-3.60</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.88</strong></td><td align="right">6.62</td><td align="right">5.47 (seed 103)</td><td align="right">-3.74</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.05</strong></td><td align="right">4.67</td><td align="right">4.25 (seed 102)</td><td align="right">-2.62</td><td align="right">36.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.01</strong></td><td align="right">1.67</td><td align="right">1.64 (seed 102)</td><td align="right">-0.66</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.84</strong></td><td align="right">2.39</td><td align="right">2.00 (seed 101)</td><td align="right">-0.55</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.12</strong></td><td align="right">2.79</td><td align="right">2.73 (seed 102)</td><td align="right">-0.67</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.50</strong></td><td align="right">3.10</td><td align="right">3.05 (seed 103)</td><td align="right">-1.60</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>2.70</strong></td><td align="right">5.67</td><td align="right">4.19 (seed 101)</td><td align="right">-2.97</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.19</strong></td><td align="right">8.96</td><td align="right">8.52 (seed 103)</td><td align="right">-3.77</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.40</strong></td><td align="right">8.22</td><td align="right">7.63 (seed 103)</td><td align="right">-3.82</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.31</strong></td><td align="right">8.40</td><td align="right">8.06 (seed 102)</td><td align="right">-3.09</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.36</strong></td><td align="right">5.69</td><td align="right">4.90 (seed 101)</td><td align="right">-2.33</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right"><strong>5.11</strong></td><td align="right">7.08</td><td align="right">6.53 (seed 103)</td><td align="right">-1.97</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>2.06</strong></td><td align="right">3.92</td><td align="right">3.29 (seed 103)</td><td align="right">-1.86</td><td align="right">4.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.14</strong></td><td align="right">5.33</td><td align="right">4.88 (seed 101)</td><td align="right">-3.19</td><td align="right">2.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>2.89</strong></td><td align="right">5.34</td><td align="right">4.90 (seed 102)</td><td align="right">-2.45</td><td align="right">4.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right"><strong>5.51</strong></td><td align="right">6.60</td><td align="right">5.98 (seed 102)</td><td align="right">-1.09</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.91</strong></td><td align="right">2.78</td><td align="right">2.05 (seed 103)</td><td align="right">-0.87</td><td align="right">8.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>2.17</strong></td><td align="right">3.02</td><td align="right">2.76 (seed 103)</td><td align="right">-0.85</td><td align="right">10.00</td><td align="right">4.67</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right"><strong>14.98</strong></td><td align="right">21.67</td><td align="right">21.34 (seed 102)</td><td align="right">-6.69</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right"><strong>21.28</strong></td><td align="right">25.13</td><td align="right">24.69 (seed 101)</td><td align="right">-3.85</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right"><strong>21.15</strong></td><td align="right">24.34</td><td align="right">23.84 (seed 102)</td><td align="right">-3.19</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.46</strong></td><td align="right"><strong>6.69</strong></td><td></td><td align="right"><strong>-2.22</strong></td><td align="right"><strong>9.78</strong></td><td align="right"><strong>7.11</strong></td><td></td></tr>
</tbody>
</table>
