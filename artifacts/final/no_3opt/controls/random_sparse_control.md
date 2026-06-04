# Random Sparse Control

This sanity check compares the final MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison reads MSA from the final results and averages the available final-control random seeds for each instance.

## Findings

- **MSA had lower average best deviation than the random-sparse mean in 22/28 instances.**
- **Mean average best deviation: MSA 3.22%, random sparse 3.78%, delta -0.56 pp.**
- **Mean success rate: MSA 5.64%, random sparse 3.71%, delta +1.93 pp.**
- MSA also beat the best random seed in 19/28 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.001514.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">4.28</td><td align="right"><strong>3.80</strong></td><td align="right">3.39 (seed 3)</td><td align="right">+0.48</td><td align="right">0.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">0.00 (seed 1)</td><td align="right">+0.00</td><td align="right">94.00</td><td align="right">85.33</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.98</strong></td><td align="right">4.39</td><td align="right">3.92 (seed 1)</td><td align="right">-0.41</td><td align="right">4.00</td><td align="right">2.00</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right">3.77</td><td align="right"><strong>3.75</strong></td><td align="right">3.58 (seed 2)</td><td align="right">+0.02</td><td align="right">0.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.63</strong></td><td align="right">3.17</td><td align="right">3.13 (seed 3)</td><td align="right">-0.54</td><td align="right">2.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.48</strong></td><td align="right">4.87</td><td align="right">4.84 (seed 3)</td><td align="right">-0.39</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_1</td><td align="right">2.97</td><td align="right"><strong>2.93</strong></td><td align="right">2.84 (seed 1)</td><td align="right">+0.04</td><td align="right">6.00</td><td align="right">3.33</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.53</strong></td><td align="right">4.89</td><td align="right">4.81 (seed 2)</td><td align="right">-2.36</td><td align="right">10.00</td><td align="right">2.00</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.05</strong></td><td align="right">2.05</td><td align="right">1.79 (seed 1)</td><td align="right">-1.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>2.05</strong></td><td align="right">3.39</td><td align="right">2.43 (seed 1)</td><td align="right">-1.34</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.18</strong></td><td align="right">2.47</td><td align="right">2.28 (seed 2)</td><td align="right">-1.29</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.39</strong></td><td align="right">3.18</td><td align="right">2.84 (seed 3)</td><td align="right">-0.79</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc188</td><td align="right"><strong>1.52</strong></td><td align="right">2.00</td><td align="right">1.74 (seed 2)</td><td align="right">-0.48</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft53</td><td align="right"><strong>5.24</strong></td><td align="right">6.48</td><td align="right">6.30 (seed 1)</td><td align="right">-1.24</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.30</strong></td><td align="right">2.81</td><td align="right">2.61 (seed 1)</td><td align="right">-0.51</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.65</strong></td><td align="right">5.22</td><td align="right">4.80 (seed 2)</td><td align="right">-1.57</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.53</strong></td><td align="right">5.98</td><td align="right">5.75 (seed 3)</td><td align="right">-0.45</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv120</td><td align="right">5.12</td><td align="right"><strong>4.62</strong></td><td align="right">4.06 (seed 3)</td><td align="right">+0.50</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.98</strong></td><td align="right">5.19</td><td align="right">5.00 (seed 2)</td><td align="right">-0.21</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.35</strong></td><td align="right">6.01</td><td align="right">5.74 (seed 1)</td><td align="right">-0.66</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv150</td><td align="right"><strong>4.26</strong></td><td align="right">4.47</td><td align="right">4.24 (seed 1)</td><td align="right">-0.21</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right">4.31</td><td align="right"><strong>4.27</strong></td><td align="right">4.10 (seed 2)</td><td align="right">+0.04</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.99</strong></td><td align="right">5.29</td><td align="right">5.05 (seed 3)</td><td align="right">-0.30</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>2.27</strong></td><td align="right">2.67</td><td align="right">2.43 (seed 2)</td><td align="right">-0.40</td><td align="right">8.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>1.01</strong></td><td align="right">1.78</td><td align="right">1.70 (seed 3)</td><td align="right">-0.77</td><td align="right">32.00</td><td align="right">5.33</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>3.17</strong></td><td align="right">3.32</td><td align="right">3.01 (seed 1)</td><td align="right">-0.15</td><td align="right">2.00</td><td align="right">2.67</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>4.33</strong></td><td align="right">5.70</td><td align="right">5.18 (seed 1)</td><td align="right">-1.37</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.82</strong></td><td align="right">1.04</td><td align="right">0.95 (seed 1)</td><td align="right">-0.22</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>3.22</strong></td><td align="right"><strong>3.78</strong></td><td></td><td align="right"><strong>-0.56</strong></td><td align="right"><strong>5.64</strong></td><td align="right"><strong>3.71</strong></td><td></td></tr>
</tbody>
</table>
