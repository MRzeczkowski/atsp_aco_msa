# MSA Count Scaling

This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.

The table uses pooled totals over 5 selected instances. Precision and recall use the 5 instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.

## Findings

- **Precision changes from 58.59% with 1 MSA to 75.13% with all MSAs.**
- **Recall changes from 30.10% with 1 MSA to 20.98% with all MSAs.**

<table>
<thead>
<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>1</td><td align="right">100.00</td><td align="right">58.59</td><td align="right">30.10</td></tr>
<tr><td>2</td><td align="right">90.14</td><td align="right">61.88</td><td align="right">28.65</td></tr>
<tr><td>4</td><td align="right">84.51</td><td align="right">63.33</td><td align="right">27.50</td></tr>
<tr><td>8</td><td align="right">76.90</td><td align="right">66.67</td><td align="right">26.34</td></tr>
<tr><td>16</td><td align="right">71.55</td><td align="right">66.54</td><td align="right">24.46</td></tr>
<tr><td>32</td><td align="right">62.25</td><td align="right">73.76</td><td align="right">23.59</td></tr>
<tr><td>64</td><td align="right">55.77</td><td align="right">75.25</td><td align="right">21.56</td></tr>
<tr><td>all</td><td align="right">54.37</td><td align="right">75.13</td><td align="right">20.98</td></tr>
</tbody>
</table>
