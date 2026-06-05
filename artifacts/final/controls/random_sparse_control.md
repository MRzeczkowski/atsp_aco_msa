# Random Sparse Control

This sanity check compares the final MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison reads MSA from the final results and averages the available final-control random seeds for each instance.

## Findings

- **MSA had lower average best deviation than the random-sparse mean in 17/26 instances.**
- **Mean average best deviation: MSA 5.57%, random sparse 5.83%, delta -0.25 pp.**
- **Mean success rate: MSA 5.38%, random sparse 4.36%, delta +1.03 pp.**
- MSA also beat the best random seed in 12/26 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.063915.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">0.00 (seed 3)</td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">98.67</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 2)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>1.94</strong></td><td align="right">2.11</td><td align="right">1.94 (seed 3)</td><td align="right">-0.17</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right">4.11</td><td align="right"><strong>3.91</strong></td><td align="right">3.30 (seed 2)</td><td align="right">+0.20</td><td align="right">4.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>3.00</strong></td><td align="right">3.25</td><td align="right">3.15 (seed 3)</td><td align="right">-0.25</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.51</strong></td><td align="right">5.07</td><td align="right">4.48 (seed 3)</td><td align="right">-0.56</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>3.50</strong></td><td align="right">4.61</td><td align="right">4.34 (seed 3)</td><td align="right">-1.11</td><td align="right">12.00</td><td align="right">3.33</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.03</strong></td><td align="right">2.12</td><td align="right">1.96 (seed 1)</td><td align="right">-1.09</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>2.21</strong></td><td align="right">3.93</td><td align="right">3.11 (seed 1)</td><td align="right">-1.72</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.48</strong></td><td align="right">3.19</td><td align="right">2.94 (seed 1)</td><td align="right">-0.71</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.19</strong></td><td align="right">2.71</td><td align="right">2.58 (seed 1)</td><td align="right">-0.52</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.50</strong></td><td align="right">4.85</td><td align="right">4.58 (seed 3)</td><td align="right">-1.35</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.59</strong></td><td align="right">5.86</td><td align="right">5.54 (seed 3)</td><td align="right">-0.27</td><td align="right">0.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.89</strong></td><td align="right">5.35</td><td align="right">5.20 (seed 2)</td><td align="right">-0.46</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.70</strong></td><td align="right">6.13</td><td align="right">5.87 (seed 1)</td><td align="right">-0.43</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.93</strong></td><td align="right">4.23</td><td align="right">3.64 (seed 1)</td><td align="right">-0.30</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right">5.40</td><td align="right"><strong>5.02</strong></td><td align="right">4.56 (seed 1)</td><td align="right">+0.38</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>1.71</strong></td><td align="right">2.28</td><td align="right">2.18 (seed 3)</td><td align="right">-0.57</td><td align="right">8.00</td><td align="right">3.33</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>1.84</strong></td><td align="right">3.77</td><td align="right">3.54 (seed 3)</td><td align="right">-1.93</td><td align="right">4.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>2.88</strong></td><td align="right">3.31</td><td align="right">2.94 (seed 2)</td><td align="right">-0.43</td><td align="right">10.00</td><td align="right">4.00</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right">5.61</td><td align="right"><strong>4.81</strong></td><td align="right">4.59 (seed 2)</td><td align="right">+0.80</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>2.86</strong></td><td align="right">3.08</td><td align="right">2.76 (seed 1)</td><td align="right">-0.22</td><td align="right">0.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right">3.58</td><td align="right"><strong>3.52</strong></td><td align="right">2.91 (seed 1)</td><td align="right">+0.06</td><td align="right">2.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right">19.68</td><td align="right"><strong>19.11</strong></td><td align="right">18.98 (seed 2)</td><td align="right">+0.57</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right">25.24</td><td align="right"><strong>23.72</strong></td><td align="right">23.22 (seed 3)</td><td align="right">+1.52</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right">27.26</td><td align="right"><strong>25.30</strong></td><td align="right">25.03 (seed 3)</td><td align="right">+1.96</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>5.57</strong></td><td align="right"><strong>5.83</strong></td><td></td><td align="right"><strong>-0.25</strong></td><td align="right"><strong>5.38</strong></td><td align="right"><strong>4.36</strong></td><td></td></tr>
</tbody>
</table>
