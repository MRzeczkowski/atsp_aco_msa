# MSA Count Scaling

This table shows how the MSA-support signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.

The table uses pooled totals over 28 selected instances. Precision and recall use the 23 instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.

## Findings

- **Precision changes from 55.21% with 1 MSA to 65.44% with all MSAs.**
- **Recall changes from 4.27% with 1 MSA to 3.18% with all MSAs.**

<table>
<thead>
<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>1</td><td align="right">100.00</td><td align="right">55.21</td><td align="right">4.27</td></tr>
<tr><td>2</td><td align="right">93.49</td><td align="right">57.66</td><td align="right">4.12</td></tr>
<tr><td>4</td><td align="right">87.56</td><td align="right">59.05</td><td align="right">3.94</td></tr>
<tr><td>8</td><td align="right">82.42</td><td align="right">60.05</td><td align="right">3.82</td></tr>
<tr><td>16</td><td align="right">79.97</td><td align="right">61.64</td><td align="right">3.66</td></tr>
<tr><td>32</td><td align="right">74.34</td><td align="right">63.28</td><td align="right">3.52</td></tr>
<tr><td>64</td><td align="right">70.24</td><td align="right">65.04</td><td align="right">3.35</td></tr>
<tr><td>all</td><td align="right">67.18</td><td align="right">65.44</td><td align="right">3.18</td></tr>
</tbody>
</table>
