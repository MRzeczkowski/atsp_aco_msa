# Shuffled MSA Control

This sanity check compares the final MSA heuristic against deterministic shuffles of the MSA mask. Each shuffle preserves the number and boost values of MSA-boosted directed edges, but assigns them to shuffled directed edges. The control uses the same `heuristicWeight=1.00`.

## Findings

- **MSA had lower average best deviation than the shuffled MSA mean in 17/26 instances.**
- **Mean average best deviation: MSA 5.57%, shuffled MSA 5.78%, delta -0.21 pp.**
- **Mean success rate: MSA 5.38%, shuffled MSA 4.62%, delta +0.77 pp.**
- MSA also beat the best shuffled MSA seed in 14/26 instances.
- Two-sided sign-test p-value for average-best-deviation wins/losses: 0.107752.

## Per-instance comparison

Negative delta means the MSA heuristic had lower average best deviation than the shuffled MSA mean.

<table>
<thead>
<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>shuffled MSA mean avg best dev. [%]</th><th>Best shuffled MSA avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>shuffled MSA success [%]</th><th>Seeds</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>0.00</strong></td><td align="right">0.00</td><td align="right">0.00 (seed 101)</td><td align="right">-0.00</td><td align="right">100.00</td><td align="right">98.00</td><td align="right">3</td></tr>
<tr><td>atex3</td><td align="right">0.27</td><td align="right"><strong>0.27</strong></td><td align="right">0.27 (seed 102)</td><td align="right">+0.00</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>atex4</td><td align="right"><strong>1.94</strong></td><td align="right">2.00</td><td align="right">1.91 (seed 103)</td><td align="right">-0.06</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane100_0</td><td align="right">4.11</td><td align="right"><strong>3.97</strong></td><td align="right">3.24 (seed 101)</td><td align="right">+0.14</td><td align="right">4.00</td><td align="right">3.33</td><td align="right">3</td></tr>
<tr><td>crane100_2</td><td align="right"><strong>3.00</strong></td><td align="right">3.54</td><td align="right">3.15 (seed 101)</td><td align="right">-0.54</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>4.51</strong></td><td align="right">5.01</td><td align="right">4.63 (seed 103)</td><td align="right">-0.50</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>3.50</strong></td><td align="right">4.88</td><td align="right">4.39 (seed 101)</td><td align="right">-1.38</td><td align="right">12.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>dc112</td><td align="right"><strong>1.03</strong></td><td align="right">2.26</td><td align="right">2.19 (seed 101)</td><td align="right">-1.23</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc126</td><td align="right"><strong>2.21</strong></td><td align="right">3.05</td><td align="right">2.81 (seed 103)</td><td align="right">-0.84</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>dc176</td><td align="right"><strong>2.48</strong></td><td align="right">3.16</td><td align="right">3.06 (seed 103)</td><td align="right">-0.68</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ft70</td><td align="right"><strong>2.19</strong></td><td align="right">2.73</td><td align="right">2.63 (seed 103)</td><td align="right">-0.54</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv100</td><td align="right"><strong>3.50</strong></td><td align="right">4.92</td><td align="right">4.14 (seed 101)</td><td align="right">-1.42</td><td align="right">0.00</td><td align="right">0.67</td><td align="right">3</td></tr>
<tr><td>ftv110</td><td align="right"><strong>5.59</strong></td><td align="right">5.97</td><td align="right">5.62 (seed 101)</td><td align="right">-0.38</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv130</td><td align="right"><strong>4.89</strong></td><td align="right">5.72</td><td align="right">5.44 (seed 103)</td><td align="right">-0.83</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv140</td><td align="right"><strong>5.70</strong></td><td align="right">6.21</td><td align="right">6.11 (seed 103)</td><td align="right">-0.51</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv160</td><td align="right"><strong>3.93</strong></td><td align="right">4.78</td><td align="right">4.55 (seed 103)</td><td align="right">-0.85</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv170</td><td align="right">5.40</td><td align="right"><strong>5.28</strong></td><td align="right">4.64 (seed 103)</td><td align="right">+0.12</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv35</td><td align="right"><strong>1.71</strong></td><td align="right">1.94</td><td align="right">1.82 (seed 102)</td><td align="right">-0.23</td><td align="right">8.00</td><td align="right">6.67</td><td align="right">3</td></tr>
<tr><td>ftv38</td><td align="right"><strong>1.84</strong></td><td align="right">3.45</td><td align="right">2.72 (seed 101)</td><td align="right">-1.61</td><td align="right">4.00</td><td align="right">2.67</td><td align="right">3</td></tr>
<tr><td>ftv44</td><td align="right"><strong>2.88</strong></td><td align="right">3.04</td><td align="right">2.82 (seed 101)</td><td align="right">-0.16</td><td align="right">10.00</td><td align="right">4.00</td><td align="right">3</td></tr>
<tr><td>ftv47</td><td align="right">5.61</td><td align="right"><strong>4.72</strong></td><td align="right">4.36 (seed 103)</td><td align="right">+0.89</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>ftv55</td><td align="right">2.86</td><td align="right"><strong>2.74</strong></td><td align="right">2.37 (seed 101)</td><td align="right">+0.12</td><td align="right">0.00</td><td align="right">2.67</td><td align="right">3</td></tr>
<tr><td>ftv70</td><td align="right">3.58</td><td align="right"><strong>3.36</strong></td><td align="right">2.79 (seed 103)</td><td align="right">+0.22</td><td align="right">2.00</td><td align="right">1.33</td><td align="right">3</td></tr>
<tr><td>rbg358</td><td align="right">19.68</td><td align="right"><strong>18.60</strong></td><td align="right">17.63 (seed 102)</td><td align="right">+1.08</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg403</td><td align="right">25.24</td><td align="right"><strong>24.08</strong></td><td align="right">23.54 (seed 101)</td><td align="right">+1.16</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td>rbg443</td><td align="right">27.26</td><td align="right"><strong>24.68</strong></td><td align="right">24.12 (seed 101)</td><td align="right">+2.58</td><td align="right">0.00</td><td align="right">0.00</td><td align="right">3</td></tr>
<tr><td><strong>Average</strong></td><td align="right"><strong>5.57</strong></td><td align="right"><strong>5.78</strong></td><td></td><td align="right"><strong>-0.21</strong></td><td align="right"><strong>5.38</strong></td><td align="right"><strong>4.62</strong></td><td></td></tr>
</tbody>
</table>
