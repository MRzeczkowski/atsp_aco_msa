# MSA vs Distance-ranked Edge Categories

This report splits directed edges into four categories: edges boosted by both the MSA-impact heuristic and the distance-ranked sparse control, edges boosted only by the MSA-impact heuristic, edges boosted only by the distance-ranked sparse control, and edges boosted by neither. The distance-ranked set uses the same boosted-edge count as the MSA-impact heuristic and the same deterministic ordering as the control heuristic.

Instances without found optimal tours in `solutions.csv` are omitted, because precision and recall need a reference edge set.

## Findings

- **Analyzed 6 instances with found optimal tours.**
- **40.56% of MSA-impact boosted edges are also distance-ranked sparse edges.**
- **MSA-only precision 62.35% vs distance-only precision 38.82%.**
- **MSA-only recall 10.67% vs distance-only recall 6.65%.**
- **Found-optimal edge distribution: both 86, MSA-only 106, distance-only 66, neither 735.**

## Pooled Categories

<table>
<thead>
<tr><th>Category</th><th>Edges</th><th>Found-optimal edges</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>Both</td><td align="right">116</td><td align="right">86</td><td align="right">74.14</td><td align="right">8.66</td></tr>
<tr><td>MSA only</td><td align="right">170</td><td align="right">106</td><td align="right">62.35</td><td align="right">10.67</td></tr>
<tr><td>Distance-ranked only</td><td align="right">170</td><td align="right">66</td><td align="right">38.82</td><td align="right">6.65</td></tr>
<tr><td>Neither</td><td align="right">39262</td><td align="right">735</td><td align="right">1.87</td><td align="right">74.02</td></tr>
</tbody>
</table>

## Per-instance Counts

<table>
<thead>
<tr><th>Instance</th><th>n</th><th>Optimal edges</th><th>Both edges</th><th>MSA-only edges</th><th>Distance-only edges</th><th>Optimal both</th><th>Optimal MSA-only</th><th>Optimal distance-only</th><th>Optimal neither</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">72</td><td align="right">360</td><td align="right">15</td><td align="right">17</td><td align="right">17</td><td align="right">15</td><td align="right">17</td><td align="right">17</td><td align="right">311</td></tr>
<tr><td>crane100_1</td><td align="right">100</td><td align="right">100</td><td align="right">36</td><td align="right">28</td><td align="right">28</td><td align="right">24</td><td align="right">11</td><td align="right">8</td><td align="right">57</td></tr>
<tr><td>ftv64</td><td align="right">65</td><td align="right">92</td><td align="right">16</td><td align="right">22</td><td align="right">22</td><td align="right">14</td><td align="right">15</td><td align="right">10</td><td align="right">53</td></tr>
<tr><td>ftv90</td><td align="right">91</td><td align="right">113</td><td align="right">21</td><td align="right">31</td><td align="right">31</td><td align="right">15</td><td align="right">23</td><td align="right">14</td><td align="right">61</td></tr>
<tr><td>ry48p</td><td align="right">48</td><td align="right">48</td><td align="right">4</td><td align="right">13</td><td align="right">13</td><td align="right">2</td><td align="right">9</td><td align="right">4</td><td align="right">33</td></tr>
<tr><td>td100_1</td><td align="right">101</td><td align="right">280</td><td align="right">24</td><td align="right">59</td><td align="right">59</td><td align="right">16</td><td align="right">31</td><td align="right">13</td><td align="right">220</td></tr>
</tbody>
</table>
