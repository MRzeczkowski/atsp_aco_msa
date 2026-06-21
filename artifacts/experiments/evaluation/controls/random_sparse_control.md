# Random Sparse Control

This sanity check compares Strict MSA against deterministic random sparse masks. Each random mask boosts the same number of directed edges as Strict MSA and uses the same heuristic weight. The comparison reads Strict MSA from the evaluation results and averages the available evaluation-control random seeds for each instance.

## Findings

- **Strict MSA had lower average best deviation than the random-sparse mean in 26/27 instances.**
- **Mean average best deviation: Strict MSA 4.49%, random sparse 6.66%, delta -2.17 pp.**
- **Mean success rate: Strict MSA 10.44%, random sparse 6.94%, delta +3.51 pp.**
- Strict MSA also beat the best random seed in 24/27 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.000000.

## Per-instance comparison

Negative delta means Strict MSA had lower average best deviation than the random-sparse mean.

<table>
<thead>
<tr><th>Instance</th><th>Strict MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>Strict MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>0.00</strong></td><td align="right">0.01</td><td align="right">0.01 (seed 2)</td><td align="right">-0.01</td><td align="right">100.00</td><td align="right">95.33</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 2)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>2.00</strong></td><td align="right">4.01</td><td align="right">3.22 (seed 3)</td><td align="right">-2.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 1)</td><td align="right">-0.00</td><td align="right">100.00</td><td align="right">82.00</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>4.07</strong></td><td align="right">6.90</td><td align="right">6.45 (seed 2)</td><td align="right">-2.83</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.76</strong></td><td align="right">6.17</td><td align="right">5.99 (seed 3)</td><td align="right">-3.41</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.40</strong></td><td align="right">5.60</td><td align="right">4.18 (seed 3)</td><td align="right">-3.20</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.04</strong></td><td align="right">5.13</td><td align="right">4.47 (seed 3)</td><td align="right">-3.09</td><td align="right">42.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.00</strong></td><td align="right">1.54</td><td align="right">1.46 (seed 1)</td><td align="right">-0.54</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.79</strong></td><td align="right">2.18</td><td align="right">2.03 (seed 3)</td><td align="right">-0.39</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.10</strong></td><td align="right">2.80</td><td align="right">2.68 (seed 3)</td><td align="right">-0.70</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.49</strong></td><td align="right">2.95</td><td align="right">2.89 (seed 2)</td><td align="right">-1.46</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.13</strong></td><td align="right">6.82</td><td align="right">5.91 (seed 2)</td><td align="right">-3.69</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.03</strong></td><td align="right">8.13</td><td align="right">7.44 (seed 3)</td><td align="right">-3.10</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.35</strong></td><td align="right">8.23</td><td align="right">7.61 (seed 3)</td><td align="right">-3.88</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>4.90</strong></td><td align="right">8.39</td><td align="right">8.10 (seed 3)</td><td align="right">-3.49</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.68</strong></td><td align="right">5.54</td><td align="right">4.86 (seed 1)</td><td align="right">-1.86</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.79</strong></td><td align="right">6.85</td><td align="right">6.56 (seed 1)</td><td align="right">-2.06</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>2.01</strong></td><td align="right">3.90</td><td align="right">3.56 (seed 3)</td><td align="right">-1.89</td><td align="right">4.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.90</strong></td><td align="right">5.16</td><td align="right">5.05 (seed 3)</td><td align="right">-2.26</td><td align="right">6.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>3.26</strong></td><td align="right">5.64</td><td align="right">4.88 (seed 2)</td><td align="right">-2.38</td><td align="right">6.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right"><strong>5.49</strong></td><td align="right">7.10</td><td align="right">6.60 (seed 2)</td><td align="right">-1.61</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.85</strong></td><td align="right">2.55</td><td align="right">2.09 (seed 1)</td><td align="right">-0.70</td><td align="right">16.00</td><td align="right">2.00</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>2.25</strong></td><td align="right">2.99</td><td align="right">2.13 (seed 1)</td><td align="right">-0.74</td><td align="right">8.00</td><td align="right">4.00</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right"><strong>15.04</strong></td><td align="right">21.59</td><td align="right">21.24 (seed 3)</td><td align="right">-6.55</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right"><strong>21.50</strong></td><td align="right">24.52</td><td align="right">24.02 (seed 3)</td><td align="right">-3.02</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right"><strong>21.15</strong></td><td align="right">24.78</td><td align="right">24.21 (seed 1)</td><td align="right">-3.63</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.49</strong></td><td align="right"><strong>6.66</strong></td><td></td><td align="right"><strong>-2.17</strong></td><td align="right"><strong>10.44</strong></td><td align="right"><strong>6.94</strong></td><td></td></tr>
</tbody>
</table>
