# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the same `heuristicWeight=0.40`.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 25/27 instances.**
- **Mean average best deviation: MSA 4.49%, shuffled MSA 6.69%, delta -2.19 pp.**
- **Mean success rate: MSA 10.44%, shuffled MSA 7.06%, delta +3.38 pp.**
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
<tr><td>atex4</td><td align="right"><strong>2.00</strong></td><td align="right">3.11</td><td align="right">2.28 (seed 102)</td><td align="right">-1.11</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>code198</td><td align="right">0.00</td><td align="right"><strong>0.00</strong></td><td align="right">0.00 (seed 103)</td><td align="right">+0.00</td><td align="right">100.00</td><td align="right">88.67</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>4.07</strong></td><td align="right">6.68</td><td align="right">6.23 (seed 103)</td><td align="right">-2.61</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>2.76</strong></td><td align="right">6.60</td><td align="right">6.24 (seed 102)</td><td align="right">-3.84</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.40</strong></td><td align="right">5.86</td><td align="right">5.26 (seed 102)</td><td align="right">-3.46</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.04</strong></td><td align="right">4.92</td><td align="right">4.03 (seed 101)</td><td align="right">-2.88</td><td align="right">42.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.00</strong></td><td align="right">1.63</td><td align="right">1.62 (seed 103)</td><td align="right">-0.63</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.79</strong></td><td align="right">2.38</td><td align="right">1.97 (seed 101)</td><td align="right">-0.59</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.10</strong></td><td align="right">2.85</td><td align="right">2.74 (seed 102)</td><td align="right">-0.75</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.49</strong></td><td align="right">2.87</td><td align="right">2.83 (seed 101)</td><td align="right">-1.38</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.13</strong></td><td align="right">6.14</td><td align="right">5.25 (seed 101)</td><td align="right">-3.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.03</strong></td><td align="right">9.12</td><td align="right">8.76 (seed 103)</td><td align="right">-4.09</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.35</strong></td><td align="right">8.01</td><td align="right">7.06 (seed 103)</td><td align="right">-3.66</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>4.90</strong></td><td align="right">8.46</td><td align="right">8.26 (seed 101)</td><td align="right">-3.56</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.68</strong></td><td align="right">5.84</td><td align="right">5.60 (seed 101)</td><td align="right">-2.16</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right"><strong>4.79</strong></td><td align="right">7.45</td><td align="right">7.41 (seed 103)</td><td align="right">-2.66</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>2.01</strong></td><td align="right">3.68</td><td align="right">3.41 (seed 103)</td><td align="right">-1.67</td><td align="right">4.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.90</strong></td><td align="right">5.40</td><td align="right">5.00 (seed 101)</td><td align="right">-2.50</td><td align="right">6.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>3.26</strong></td><td align="right">4.84</td><td align="right">4.24 (seed 102)</td><td align="right">-1.58</td><td align="right">6.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right"><strong>5.49</strong></td><td align="right">7.37</td><td align="right">7.03 (seed 101)</td><td align="right">-1.88</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.85</strong></td><td align="right">2.80</td><td align="right">2.17 (seed 103)</td><td align="right">-0.95</td><td align="right">16.00</td><td align="right">2.00</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right"><strong>2.25</strong></td><td align="right">3.03</td><td align="right">2.88 (seed 103)</td><td align="right">-0.78</td><td align="right">8.00</td><td align="right">2.67</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right"><strong>15.04</strong></td><td align="right">21.75</td><td align="right">21.47 (seed 102)</td><td align="right">-6.71</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right"><strong>21.50</strong></td><td align="right">25.10</td><td align="right">24.60 (seed 101)</td><td align="right">-3.60</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right"><strong>21.15</strong></td><td align="right">24.33</td><td align="right">23.84 (seed 101)</td><td align="right">-3.18</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.49</strong></td><td align="right"><strong>6.69</strong></td><td></td><td align="right"><strong>-2.19</strong></td><td align="right"><strong>10.44</strong></td><td align="right"><strong>7.06</strong></td><td></td></tr>
</tbody>
</table>
