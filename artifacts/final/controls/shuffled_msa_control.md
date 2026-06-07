# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the same `heuristicWeight=0.70`.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 16/26 instances.**
- **Mean average best deviation: MSA 4.60%, shuffled MSA 4.79%, delta -0.19 pp.**
- **Mean success rate: MSA 7.69%, shuffled MSA 5.92%, delta +1.77 pp.**
- MSA also beat the best shuffled MSA seed in 13/26 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.229523.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>shuffled MSA mean avg best dev. [%]</th><th>Best shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 103)</td><td align="right">-0.00</td><td align="right">98.00</td><td align="right">98.00</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 102)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>1.97</strong></td><td align="right">1.98</td><td align="right">1.88 (seed 103)</td><td align="right">-0.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right"><strong>3.62</strong></td><td align="right">4.25</td><td align="right">4.03 (seed 101)</td><td align="right">-0.63</td><td align="right">2.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right">3.09</td><td align="right"><strong>2.96</strong></td><td align="right">2.59 (seed 102)</td><td align="right">+0.13</td><td align="right">2.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>2.28</strong></td><td align="right">2.97</td><td align="right">2.83 (seed 103)</td><td align="right">-0.69</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>2.05</strong></td><td align="right">3.04</td><td align="right">2.86 (seed 102)</td><td align="right">-0.99</td><td align="right">42.00</td><td align="right">22.00</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.02</strong></td><td align="right">2.08</td><td align="right">1.83 (seed 101)</td><td align="right">-1.06</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>1.94</strong></td><td align="right">2.59</td><td align="right">2.28 (seed 101)</td><td align="right">-0.65</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.28</strong></td><td align="right">2.53</td><td align="right">2.37 (seed 103)</td><td align="right">-0.25</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>1.56</strong></td><td align="right">1.76</td><td align="right">1.69 (seed 103)</td><td align="right">-0.20</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.40</strong></td><td align="right">5.06</td><td align="right">4.88 (seed 101)</td><td align="right">-1.66</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>4.42</strong></td><td align="right">6.43</td><td align="right">6.25 (seed 102)</td><td align="right">-2.01</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.44</strong></td><td align="right">4.78</td><td align="right">4.48 (seed 103)</td><td align="right">-0.34</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>4.85</strong></td><td align="right">4.93</td><td align="right">4.79 (seed 103)</td><td align="right">-0.08</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right">3.15</td><td align="right"><strong>3.00</strong></td><td align="right">2.63 (seed 102)</td><td align="right">+0.15</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right">4.83</td><td align="right"><strong>4.19</strong></td><td align="right">4.02 (seed 102)</td><td align="right">+0.64</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right">2.20</td><td align="right"><strong>2.01</strong></td><td align="right">1.90 (seed 103)</td><td align="right">+0.19</td><td align="right">4.00</td><td align="right">6.67</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>2.36</strong></td><td align="right">3.43</td><td align="right">3.34 (seed 102)</td><td align="right">-1.07</td><td align="right">2.00</td><td align="right">2.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>2.88</strong></td><td align="right">3.33</td><td align="right">2.96 (seed 101)</td><td align="right">-0.45</td><td align="right">12.00</td><td align="right">5.33</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right">5.58</td><td align="right"><strong>4.72</strong></td><td align="right">4.41 (seed 103)</td><td align="right">+0.86</td><td align="right">0.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right"><strong>1.50</strong></td><td align="right">2.02</td><td align="right">1.88 (seed 103)</td><td align="right">-0.52</td><td align="right">20.00</td><td align="right">6.67</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right">2.40</td><td align="right"><strong>2.39</strong></td><td align="right">1.86 (seed 103)</td><td align="right">+0.01</td><td align="right">18.00</td><td align="right">9.33</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right">14.84</td><td align="right"><strong>14.15</strong></td><td align="right">13.83 (seed 102)</td><td align="right">+0.69</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right">21.51</td><td align="right"><strong>20.44</strong></td><td align="right">20.32 (seed 101)</td><td align="right">+1.07</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right">21.07</td><td align="right"><strong>19.12</strong></td><td align="right">18.87 (seed 102)</td><td align="right">+1.95</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>4.60</strong></td><td align="right"><strong>4.79</strong></td><td></td><td align="right"><strong>-0.19</strong></td><td align="right"><strong>7.69</strong></td><td align="right"><strong>5.92</strong></td><td></td></tr>
</tbody>
</table>
