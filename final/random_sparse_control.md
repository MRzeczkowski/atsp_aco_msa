# Random Sparse Control

This sanity check compares the real MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison uses MSA `heuristicWeight=0.90` and averages the available random seeds for each instance.

## Findings

- **MSA had lower average best deviation than the random-sparse mean in 24/28 instances.**
- **Mean average best deviation: MSA 3.21%, random sparse 3.72%, delta -0.51 pp.**
- **Mean success rate: MSA 5.71%, random sparse 3.57%, delta +2.14 pp.**
- MSA also beat the best random seed in 17/28 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.000180.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">4.37</td><td align="right"><strong>3.85</strong></td><td align="right">3.31 (seed 3)</td><td align="right">+0.52</td><td align="right">0.00</td><td align="right">1.11</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 1)</td><td align="right">-0.00</td><td align="right">93.33</td><td align="right">82.22</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>4.52</strong></td><td align="right">4.56</td><td align="right">3.94 (seed 1)</td><td align="right">-0.04</td><td align="right">0.00</td><td align="right">2.22</td><td align="right">3</td></tr>
<tr><td>crane100_1</td><td align="right"><strong>3.74</strong></td><td align="right">3.76</td><td align="right">3.69 (seed 1)</td><td align="right">-0.02</td><td align="right">0.00</td><td align="right">2.22</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.70</strong></td><td align="right">3.18</td><td align="right">3.03 (seed 1)</td><td align="right">-0.48</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.44</strong></td><td align="right">4.98</td><td align="right">4.65 (seed 1)</td><td align="right">-0.54</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_1</td><td align="right"><strong>3.10</strong></td><td align="right">3.14</td><td align="right">3.06 (seed 1)</td><td align="right">-0.04</td><td align="right">3.33</td><td align="right">2.22</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.60</strong></td><td align="right">4.74</td><td align="right">4.62 (seed 2)</td><td align="right">-2.14</td><td align="right">13.33</td><td align="right">1.11</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.05</strong></td><td align="right">2.06</td><td align="right">1.76 (seed 1)</td><td align="right">-1.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>2.08</strong></td><td align="right">3.40</td><td align="right">2.41 (seed 1)</td><td align="right">-1.32</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc134</td><td align="right"><strong>1.17</strong></td><td align="right">2.51</td><td align="right">2.32 (seed 2)</td><td align="right">-1.34</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.39</strong></td><td align="right">3.17</td><td align="right">2.85 (seed 3)</td><td align="right">-0.78</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc188</td><td align="right"><strong>1.51</strong></td><td align="right">1.99</td><td align="right">1.72 (seed 2)</td><td align="right">-0.48</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft53</td><td align="right"><strong>4.90</strong></td><td align="right">6.46</td><td align="right">5.75 (seed 2)</td><td align="right">-1.56</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.24</strong></td><td align="right">2.71</td><td align="right">2.47 (seed 1)</td><td align="right">-0.47</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.62</strong></td><td align="right">4.81</td><td align="right">4.38 (seed 2)</td><td align="right">-1.19</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.72</strong></td><td align="right">5.84</td><td align="right">5.46 (seed 2)</td><td align="right">-0.12</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv120</td><td align="right">4.90</td><td align="right"><strong>4.58</strong></td><td align="right">3.99 (seed 3)</td><td align="right">+0.32</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>5.12</strong></td><td align="right">5.13</td><td align="right">4.99 (seed 2)</td><td align="right">-0.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.03</strong></td><td align="right">5.83</td><td align="right">5.48 (seed 1)</td><td align="right">-0.80</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv150</td><td align="right">4.52</td><td align="right"><strong>4.42</strong></td><td align="right">4.21 (seed 3)</td><td align="right">+0.10</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right">4.09</td><td align="right"><strong>4.00</strong></td><td align="right">3.64 (seed 3)</td><td align="right">+0.09</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.94</strong></td><td align="right">5.14</td><td align="right">4.93 (seed 3)</td><td align="right">-0.20</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>2.32</strong></td><td align="right">2.53</td><td align="right">2.41 (seed 2)</td><td align="right">-0.21</td><td align="right">10.00</td><td align="right">1.11</td><td align="right">3</td></tr>
<tr><td>ftv64</td><td align="right"><strong>1.07</strong></td><td align="right">1.66</td><td align="right">1.48 (seed 3)</td><td align="right">-0.59</td><td align="right">36.67</td><td align="right">6.67</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>2.92</strong></td><td align="right">3.31</td><td align="right">2.97 (seed 1)</td><td align="right">-0.39</td><td align="right">3.33</td><td align="right">1.11</td><td align="right">3</td></tr>
<tr><td>ftv90</td><td align="right"><strong>4.04</strong></td><td align="right">5.39</td><td align="right">4.91 (seed 1)</td><td align="right">-1.35</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>td100_1</td><td align="right"><strong>0.86</strong></td><td align="right">1.05</td><td align="right">1.00 (seed 2)</td><td align="right">-0.19</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>3.21</strong></td><td align="right"><strong>3.72</strong></td><td></td><td align="right"><strong>-0.51</strong></td><td align="right"><strong>5.71</strong></td><td align="right"><strong>3.57</strong></td><td></td></tr>
</tbody>
</table>
