# MSA Count Scaling

This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.

The table uses pooled totals over 7 selected instances. Precision and recall use the 6 instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.

## Findings

- **Precision changes from 56.69% with 1 MSA to 67.13% with all MSAs.**
- **Recall changes from 26.89% with 1 MSA to 19.34% with all MSAs.**

<table>
<thead>
<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>1</td><td align="right">100.00</td><td align="right">56.69</td><td align="right">26.89</td></tr>
<tr><td>2</td><td align="right">92.22</td><td align="right">60.38</td><td align="right">25.78</td></tr>
<tr><td>4</td><td align="right">89.40</td><td align="right">60.29</td><td align="right">24.77</td></tr>
<tr><td>8</td><td align="right">85.10</td><td align="right">62.57</td><td align="right">24.07</td></tr>
<tr><td>16</td><td align="right">80.46</td><td align="right">63.20</td><td align="right">22.66</td></tr>
<tr><td>32</td><td align="right">75.17</td><td align="right">66.56</td><td align="right">21.65</td></tr>
<tr><td>64</td><td align="right">70.03</td><td align="right">67.24</td><td align="right">19.84</td></tr>
<tr><td>all</td><td align="right">68.71</td><td align="right">67.13</td><td align="right">19.34</td></tr>
</tbody>
</table>
