# Random Sparse Control

This sanity check compares the final MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison reads MSA from the final results and averages the available final-control random seeds for each instance.

## Findings

- **MSA had lower average best deviation than the random-sparse mean in 17/26 instances.**
- **Mean average best deviation: MSA 4.60%, random sparse 4.83%, delta -0.23 pp.**
- **Mean success rate: MSA 7.69%, random sparse 5.21%, delta +2.49 pp.**
- MSA also beat the best random seed in 12/26 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.063915.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">0.00 (seed 1)</td><td align="right">+0.00</td><td align="right">98.00</td><td align="right">98.67</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 2)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>1.97</strong></td><td align="right">1.99</td><td align="right">1.92 (seed 3)</td><td align="right">-0.02</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.62</strong></td><td align="right">4.06</td><td align="right">3.61 (seed 2)</td><td align="right">-0.44</td><td align="right">2.00</td><td align="right">2.67</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right">3.09</td><td align="right"><strong>2.88</strong></td><td align="right">2.28 (seed 3)</td><td align="right">+0.21</td><td align="right">2.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.28</strong></td><td align="right">2.90</td><td align="right">2.40 (seed 3)</td><td align="right">-0.62</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.05</strong></td><td align="right">3.16</td><td align="right">2.74 (seed 3)</td><td align="right">-1.11</td><td align="right">42.00</td><td align="right">21.33</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.02</strong></td><td align="right">1.87</td><td align="right">1.70 (seed 1)</td><td align="right">-0.85</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.94</strong></td><td align="right">2.79</td><td align="right">2.07 (seed 1)</td><td align="right">-0.85</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.28</strong></td><td align="right">2.71</td><td align="right">2.49 (seed 3)</td><td align="right">-0.43</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.56</strong></td><td align="right">1.72</td><td align="right">1.65 (seed 1)</td><td align="right">-0.16</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.40</strong></td><td align="right">5.32</td><td align="right">5.28 (seed 2)</td><td align="right">-1.92</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>4.42</strong></td><td align="right">5.52</td><td align="right">5.47 (seed 3)</td><td align="right">-1.10</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.44</strong></td><td align="right">4.78</td><td align="right">4.64 (seed 2)</td><td align="right">-0.34</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>4.85</strong></td><td align="right">4.99</td><td align="right">4.77 (seed 2)</td><td align="right">-0.14</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right">3.15</td><td align="right"><strong>3.11</strong></td><td align="right">2.81 (seed 2)</td><td align="right">+0.04</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right">4.83</td><td align="right"><strong>4.36</strong></td><td align="right">4.21 (seed 1)</td><td align="right">+0.47</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>2.20</strong></td><td align="right">2.50</td><td align="right">2.16 (seed 3)</td><td align="right">-0.30</td><td align="right">4.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.36</strong></td><td align="right">3.57</td><td align="right">3.43 (seed 3)</td><td align="right">-1.21</td><td align="right">2.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>2.88</strong></td><td align="right">3.64</td><td align="right">3.42 (seed 2)</td><td align="right">-0.76</td><td align="right">12.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right">5.58</td><td align="right"><strong>4.84</strong></td><td align="right">4.60 (seed 2)</td><td align="right">+0.74</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.50</strong></td><td align="right">2.20</td><td align="right">2.00 (seed 1)</td><td align="right">-0.70</td><td align="right">20.00</td><td align="right">4.00</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>2.40</strong></td><td align="right">2.55</td><td align="right">2.16 (seed 2)</td><td align="right">-0.15</td><td align="right">18.00</td><td align="right">4.67</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right">14.84</td><td align="right"><strong>14.31</strong></td><td align="right">14.13 (seed 3)</td><td align="right">+0.53</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right">21.51</td><td align="right"><strong>19.93</strong></td><td align="right">19.20 (seed 3)</td><td align="right">+1.58</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right">21.07</td><td align="right"><strong>19.53</strong></td><td align="right">19.12 (seed 3)</td><td align="right">+1.54</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.60</strong></td><td align="right"><strong>4.83</strong></td><td></td><td align="right"><strong>-0.23</strong></td><td align="right"><strong>7.69</strong></td><td align="right"><strong>5.21</strong></td><td></td></tr>
</tbody>
</table>
