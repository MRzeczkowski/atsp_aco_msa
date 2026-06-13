# MSA vs Distance-ranked Edge Categories

This report splits directed edges into four categories: edges boosted by both the MSA-impact heuristic and the distance-ranked sparse control, edges boosted only by the MSA-impact heuristic, edges boosted only by the distance-ranked sparse control, and edges boosted by neither. The distance-ranked set uses the same boosted-edge count as the MSA-impact heuristic and the same deterministic ordering as the control heuristic.

Instances without found optimal tours in `solutions.csv` are omitted, because precision and recall need a reference edge set.

## Findings

- **Analyzed 6 instances with found optimal tours.**
- **48.62% of MSA-impact boosted edges are also distance-ranked sparse edges.**
- **MSA-only precision 56.20% vs distance-only precision 39.26%.**
- **MSA-only recall 13.70% vs distance-only recall 9.57%.**
- **Found-optimal edge distribution: both 156, MSA-only 136, distance-only 95, neither 606.**

## Pooled Categories

<table>
<thead>
<tr><th>Category</th><th>Edges</th><th>Found-optimal edges</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>Both</td><td align="right">229</td><td align="right">156</td><td align="right">68.12</td><td align="right">15.71</td></tr>
<tr><td>MSA only</td><td align="right">242</td><td align="right">136</td><td align="right">56.20</td><td align="right">13.70</td></tr>
<tr><td>Distance-ranked only</td><td align="right">242</td><td align="right">95</td><td align="right">39.26</td><td align="right">9.57</td></tr>
<tr><td>Neither</td><td align="right">39005</td><td align="right">606</td><td align="right">1.55</td><td align="right">61.03</td></tr>
</tbody>
</table>

## Per-instance Counts

<table>
<thead>
<tr><th>Instance</th><th>n</th><th>Optimal edges</th><th>Both edges</th><th>MSA-only edges</th><th>Distance-only edges</th><th>Optimal both</th><th>Optimal MSA-only</th><th>Optimal distance-only</th><th>Optimal neither</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">72</td><td align="right">360</td><td align="right">31</td><td align="right">40</td><td align="right">40</td><td align="right">31</td><td align="right">32</td><td align="right">40</td><td align="right">257</td></tr>
<tr><td>crane100_1</td><td align="right">100</td><td align="right">100</td><td align="right">56</td><td align="right">43</td><td align="right">43</td><td align="right">35</td><td align="right">14</td><td align="right">6</td><td align="right">45</td></tr>
<tr><td>ftv64</td><td align="right">65</td><td align="right">92</td><td align="right">38</td><td align="right">26</td><td align="right">26</td><td align="right">25</td><td align="right">15</td><td align="right">11</td><td align="right">41</td></tr>
<tr><td>ftv90</td><td align="right">91</td><td align="right">113</td><td align="right">52</td><td align="right">38</td><td align="right">38</td><td align="right">34</td><td align="right">23</td><td align="right">16</td><td align="right">40</td></tr>
<tr><td>ry48p</td><td align="right">48</td><td align="right">48</td><td align="right">19</td><td align="right">28</td><td align="right">28</td><td align="right">12</td><td align="right">17</td><td align="right">6</td><td align="right">13</td></tr>
<tr><td>td100_1</td><td align="right">101</td><td align="right">280</td><td align="right">33</td><td align="right">67</td><td align="right">67</td><td align="right">19</td><td align="right">35</td><td align="right">16</td><td align="right">210</td></tr>
</tbody>
</table>
