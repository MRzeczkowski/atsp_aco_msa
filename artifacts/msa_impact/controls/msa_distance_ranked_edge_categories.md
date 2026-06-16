# MSA vs Distance-ranked Edge Categories

This report splits directed edges into four categories: edges boosted by at least one rooted MSA-impact matrix and at least one rooted distance-ranked sparse control matrix, edges boosted only by rooted MSA, edges boosted only by rooted distance-ranked sparse control, and edges boosted by neither. The distance-ranked set uses the same per-root boosted-edge count as the rooted MSA heuristic and the same deterministic ordering as the control heuristic.

Instances without found optimal tours in `solutions.csv` are omitted, because precision and recall need a reference edge set.

## Findings

- **Analyzed 6 instances with found optimal tours.**
- **40.66% of MSA-impact boosted edges are also distance-ranked sparse edges.**
- **MSA-only precision 49.48% vs distance-only precision 34.29%.**
- **MSA-only recall 24.04% vs distance-only recall 4.83%.**
- **Found-optimal edge distribution: both 203, MSA-only 239, distance-only 48, neither 504.**

## Pooled Categories

<table>
<thead>
<tr><th>Category</th><th>Edges</th><th>Found-optimal edges</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>Both</td><td align="right">331</td><td align="right">203</td><td align="right">61.33</td><td align="right">20.42</td></tr>
<tr><td>MSA only</td><td align="right">483</td><td align="right">239</td><td align="right">49.48</td><td align="right">24.04</td></tr>
<tr><td>Distance-ranked only</td><td align="right">140</td><td align="right">48</td><td align="right">34.29</td><td align="right">4.83</td></tr>
<tr><td>Neither</td><td align="right">38764</td><td align="right">504</td><td align="right">1.30</td><td align="right">50.70</td></tr>
</tbody>
</table>

## Per-instance Counts

<table>
<thead>
<tr><th>Instance</th><th>n</th><th>Optimal edges</th><th>Both edges</th><th>MSA-only edges</th><th>Distance-only edges</th><th>Optimal both</th><th>Optimal MSA-only</th><th>Optimal distance-only</th><th>Optimal neither</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right">72</td><td align="right">361</td><td align="right">49</td><td align="right">91</td><td align="right">22</td><td align="right">49</td><td align="right">77</td><td align="right">22</td><td align="right">213</td></tr>
<tr><td>crane100_1</td><td align="right">100</td><td align="right">100</td><td align="right">70</td><td align="right">97</td><td align="right">29</td><td align="right">36</td><td align="right">29</td><td align="right">5</td><td align="right">30</td></tr>
<tr><td>ftv64</td><td align="right">65</td><td align="right">92</td><td align="right">52</td><td align="right">47</td><td align="right">12</td><td align="right">35</td><td align="right">24</td><td align="right">1</td><td align="right">32</td></tr>
<tr><td>ftv90</td><td align="right">91</td><td align="right">113</td><td align="right">81</td><td align="right">60</td><td align="right">9</td><td align="right">45</td><td align="right">32</td><td align="right">5</td><td align="right">31</td></tr>
<tr><td>ry48p</td><td align="right">48</td><td align="right">48</td><td align="right">32</td><td align="right">58</td><td align="right">15</td><td align="right">16</td><td align="right">21</td><td align="right">2</td><td align="right">9</td></tr>
<tr><td>td100_1</td><td align="right">101</td><td align="right">280</td><td align="right">47</td><td align="right">130</td><td align="right">53</td><td align="right">22</td><td align="right">56</td><td align="right">13</td><td align="right">189</td></tr>
</tbody>
</table>
