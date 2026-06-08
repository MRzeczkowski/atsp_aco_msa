# MSA vs Distance-ranked Edge Categories

This report splits directed edges into four categories: edges boosted by both strict MSA and the distance-ranked sparse control, edges boosted only by strict MSA, edges boosted only by the distance-ranked sparse control, and edges boosted by neither. The distance-ranked set uses the same boosted-edge count as strict MSA and the same deterministic ordering as the control heuristic.

Instances without found optimal tours in `solutions.csv` are omitted, because precision and recall need a reference edge set.

## Findings

- **Analyzed 7 instances with found optimal tours.**
- **61.51% of strict MSA boosted edges are also distance-ranked sparse edges.**
- **MSA-only precision 63.13% vs distance-only precision 37.43%.**
- **MSA-only recall 0.40% vs distance-only recall 0.23%.**
- **Found-optimal edge distribution: both 213, MSA-only 113, distance-only 67, neither 28141.**

## Pooled Categories

<table>
<thead>
<tr><th>Category</th><th>Edges</th><th>Found-optimal edges</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>Both</td><td align="right">286</td><td align="right">213</td><td align="right">74.48</td><td align="right">0.75</td></tr>
<tr><td>MSA only</td><td align="right">179</td><td align="right">113</td><td align="right">63.13</td><td align="right">0.40</td></tr>
<tr><td>Distance-ranked only</td><td align="right">179</td><td align="right">67</td><td align="right">37.43</td><td align="right">0.23</td></tr>
<tr><td>Neither</td><td align="right">78080</td><td align="right">28141</td><td align="right">36.04</td><td align="right">98.62</td></tr>
</tbody>
</table>

## Per-instance Counts

<table>
<thead>
<tr><th>Instance</th><th>n</th><th>Optimal edges</th><th>Both edges</th><th>MSA-only edges</th><th>Distance-only edges</th><th>Optimal both</th><th>Optimal MSA-only</th><th>Optimal distance-only</th><th>Optimal neither</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">72</td><td align="right">360</td><td align="right">15</td><td align="right">17</td><td align="right">17</td><td align="right">15</td><td align="right">17</td><td align="right">17</td><td align="right">311</td></tr>
<tr><td>code198</td><td align="right">198</td><td align="right">27541</td><td align="right">176</td><td align="right">11</td><td align="right">11</td><td align="right">130</td><td align="right">6</td><td align="right">2</td><td align="right">27403</td></tr>
<tr><td>crane100_1</td><td align="right">100</td><td align="right">100</td><td align="right">31</td><td align="right">26</td><td align="right">26</td><td align="right">21</td><td align="right">12</td><td align="right">8</td><td align="right">59</td></tr>
<tr><td>ftv64</td><td align="right">65</td><td align="right">92</td><td align="right">16</td><td align="right">22</td><td align="right">22</td><td align="right">14</td><td align="right">15</td><td align="right">10</td><td align="right">53</td></tr>
<tr><td>ftv90</td><td align="right">91</td><td align="right">113</td><td align="right">21</td><td align="right">31</td><td align="right">31</td><td align="right">15</td><td align="right">23</td><td align="right">14</td><td align="right">61</td></tr>
<tr><td>ry48p</td><td align="right">48</td><td align="right">48</td><td align="right">4</td><td align="right">13</td><td align="right">13</td><td align="right">2</td><td align="right">9</td><td align="right">4</td><td align="right">33</td></tr>
<tr><td>td100_1</td><td align="right">101</td><td align="right">280</td><td align="right">23</td><td align="right">59</td><td align="right">59</td><td align="right">16</td><td align="right">31</td><td align="right">12</td><td align="right">221</td></tr>
</tbody>
</table>
