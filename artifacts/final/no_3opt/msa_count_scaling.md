# MSA Count Scaling

This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.

The table uses pooled totals over 26 selected instances. Precision and recall use the 23 instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.

## Findings

- **Precision changes from 59.35% with 1 MSA to 64.11% with all MSAs.**
- **Recall changes from 8.07% with 1 MSA to 6.10% with all MSAs.**

<table>
<thead>
<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>1</td><td align="right">100.00</td><td align="right">59.35</td><td align="right">8.07</td></tr>
<tr><td>2</td><td align="right">95.38</td><td align="right">60.75</td><td align="right">7.87</td></tr>
<tr><td>4</td><td align="right">89.44</td><td align="right">61.99</td><td align="right">7.57</td></tr>
<tr><td>8</td><td align="right">87.22</td><td align="right">62.60</td><td align="right">7.45</td></tr>
<tr><td>16</td><td align="right">82.97</td><td align="right">63.79</td><td align="right">7.17</td></tr>
<tr><td>32</td><td align="right">80.01</td><td align="right">64.25</td><td align="right">6.94</td></tr>
<tr><td>64</td><td align="right">76.82</td><td align="right">64.97</td><td align="right">6.70</td></tr>
<tr><td>all</td><td align="right">71.54</td><td align="right">64.11</td><td align="right">6.10</td></tr>
</tbody>
</table>
