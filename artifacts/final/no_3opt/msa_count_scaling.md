# MSA Count Scaling

This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.

The table uses pooled totals over 16 selected instances. Precision and recall use the 13 instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.

## Findings

- **Precision changes from 53.36% with 1 MSA to 61.39% with all MSAs.**
- **Recall changes from 2.09% with 1 MSA to 1.56% with all MSAs.**

<table>
<thead>
<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>1</td><td align="right">99.89</td><td align="right">53.36</td><td align="right">2.09</td></tr>
<tr><td>2</td><td align="right">94.70</td><td align="right">55.57</td><td align="right">2.00</td></tr>
<tr><td>4</td><td align="right">91.87</td><td align="right">56.14</td><td align="right">1.93</td></tr>
<tr><td>8</td><td align="right">77.26</td><td align="right">57.40</td><td align="right">1.85</td></tr>
<tr><td>16</td><td align="right">77.48</td><td align="right">58.16</td><td align="right">1.78</td></tr>
<tr><td>32</td><td align="right">70.26</td><td align="right">60.22</td><td align="right">1.70</td></tr>
<tr><td>64</td><td align="right">66.53</td><td align="right">61.09</td><td align="right">1.59</td></tr>
<tr><td>all</td><td align="right">64.11</td><td align="right">61.39</td><td align="right">1.56</td></tr>
</tbody>
</table>
