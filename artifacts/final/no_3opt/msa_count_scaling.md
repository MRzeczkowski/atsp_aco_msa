# MSA Count Scaling

This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.

The table uses pooled totals over 27 selected instances. Precision and recall use the 24 instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.

## Findings

- **Precision changes from 60.03% with 1 MSA to 65.07% with all MSAs.**
- **Recall changes from 3.73% with 1 MSA to 2.90% with all MSAs.**

<table>
<thead>
<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>1</td><td align="right">100.00</td><td align="right">60.03</td><td align="right">3.73</td></tr>
<tr><td>2</td><td align="right">96.85</td><td align="right">61.46</td><td align="right">3.69</td></tr>
<tr><td>4</td><td align="right">91.73</td><td align="right">62.81</td><td align="right">3.57</td></tr>
<tr><td>8</td><td align="right">88.99</td><td align="right">63.34</td><td align="right">3.48</td></tr>
<tr><td>16</td><td align="right">85.31</td><td align="right">64.52</td><td align="right">3.38</td></tr>
<tr><td>32</td><td align="right">82.22</td><td align="right">64.82</td><td align="right">3.26</td></tr>
<tr><td>64</td><td align="right">78.99</td><td align="right">65.76</td><td align="right">3.16</td></tr>
<tr><td>all</td><td align="right">73.96</td><td align="right">65.07</td><td align="right">2.90</td></tr>
</tbody>
</table>
