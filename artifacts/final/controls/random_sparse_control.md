# Random Sparse Control

This sanity check compares the final MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison reads MSA from the final results and averages the available final-control random seeds for each instance.

## Findings

- **MSA had lower average best deviation than the random-sparse mean in 26/27 instances.**
- **Mean average best deviation: MSA 4.46%, random sparse 6.66%, delta -2.19 pp.**
- **Mean success rate: MSA 9.78%, random sparse 6.94%, delta +2.84 pp.**
- MSA also beat the best random seed in 24/27 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.000000.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>0.00</strong></td><td align="right">0.01</td><td align="right">0.01 (seed 2)</td><td align="right">-0.01</td><td align="right">100.00</td><td align="right">95.33</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 2)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>2.00</strong></td><td align="right">4.01</td><td align="right">3.22 (seed 3)</td><td align="right">-2.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 1)</td><td align="right">-0.00</td><td align="right">100.00</td><td align="right">82.00</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.91</strong></td><td align="right">6.90</td><td align="right">6.45 (seed 2)</td><td align="right">-2.99</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.79</strong></td><td align="right">6.17</td><td align="right">5.99 (seed 3)</td><td align="right">-3.38</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.88</strong></td><td align="right">5.60</td><td align="right">4.18 (seed 3)</td><td align="right">-2.72</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.05</strong></td><td align="right">5.13</td><td align="right">4.47 (seed 3)</td><td align="right">-3.08</td><td align="right">36.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.01</strong></td><td align="right">1.54</td><td align="right">1.46 (seed 1)</td><td align="right">-0.53</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.84</strong></td><td align="right">2.18</td><td align="right">2.03 (seed 3)</td><td align="right">-0.34</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.12</strong></td><td align="right">2.80</td><td align="right">2.68 (seed 3)</td><td align="right">-0.68</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.50</strong></td><td align="right">2.95</td><td align="right">2.89 (seed 2)</td><td align="right">-1.45</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>2.70</strong></td><td align="right">6.82</td><td align="right">5.91 (seed 2)</td><td align="right">-4.12</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.19</strong></td><td align="right">8.13</td><td align="right">7.44 (seed 3)</td><td align="right">-2.94</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.40</strong></td><td align="right">8.23</td><td align="right">7.61 (seed 3)</td><td align="right">-3.83</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.31</strong></td><td align="right">8.39</td><td align="right">8.10 (seed 3)</td><td align="right">-3.08</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.36</strong></td><td align="right">5.54</td><td align="right">4.86 (seed 1)</td><td align="right">-2.18</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right"><strong>5.11</strong></td><td align="right">6.85</td><td align="right">6.56 (seed 1)</td><td align="right">-1.74</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>2.06</strong></td><td align="right">3.90</td><td align="right">3.56 (seed 3)</td><td align="right">-1.84</td><td align="right">4.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.14</strong></td><td align="right">5.16</td><td align="right">5.05 (seed 3)</td><td align="right">-3.02</td><td align="right">2.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>2.89</strong></td><td align="right">5.64</td><td align="right">4.88 (seed 2)</td><td align="right">-2.75</td><td align="right">4.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right"><strong>5.51</strong></td><td align="right">7.10</td><td align="right">6.60 (seed 2)</td><td align="right">-1.59</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.91</strong></td><td align="right">2.55</td><td align="right">2.09 (seed 1)</td><td align="right">-0.64</td><td align="right">8.00</td><td align="right">2.00</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>2.17</strong></td><td align="right">2.99</td><td align="right">2.13 (seed 1)</td><td align="right">-0.82</td><td align="right">10.00</td><td align="right">4.00</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right"><strong>14.98</strong></td><td align="right">21.59</td><td align="right">21.24 (seed 3)</td><td align="right">-6.61</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right"><strong>21.28</strong></td><td align="right">24.52</td><td align="right">24.02 (seed 3)</td><td align="right">-3.24</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right"><strong>21.15</strong></td><td align="right">24.78</td><td align="right">24.21 (seed 1)</td><td align="right">-3.63</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.46</strong></td><td align="right"><strong>6.66</strong></td><td></td><td align="right"><strong>-2.19</strong></td><td align="right"><strong>9.78</strong></td><td align="right"><strong>6.94</strong></td><td></td></tr>
</tbody>
</table>
