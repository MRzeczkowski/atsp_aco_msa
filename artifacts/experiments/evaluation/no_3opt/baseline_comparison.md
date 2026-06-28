# Baseline Comparison Summary

This report compares every heuristic against the baseline MMAS variant. Negative average-best-deviation delta means the heuristic had lower deviation than the baseline.

The `rbg` split is included because this family can strongly affect aggregate averages.

## All

Instances: 27

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. delta [pp]</th><th>Wins</th><th>Ties</th><th>Losses</th><th>Sign-test p-value</th><th>Success delta [pp]</th></tr>
</thead>
<tbody>
<tr><td>Strict MSA</td><td align="right">-0.11</td><td align="right">13</td><td align="right">3</td><td align="right">11</td><td align="right">0.838820</td><td align="right">+0.89</td></tr>
<tr><td>Rooted MSA</td><td align="right">+0.05</td><td align="right">13</td><td align="right">3</td><td align="right">11</td><td align="right">0.838820</td><td align="right">+1.26</td></tr>
<tr><td>Cycle cover</td><td align="right">-1.81</td><td align="right">17</td><td align="right">3</td><td align="right">7</td><td align="right">0.063915</td><td align="right">-2.15</td></tr>
<tr><td>Cycle-cover patching</td><td align="right">-2.66</td><td align="right">19</td><td align="right">1</td><td align="right">7</td><td align="right">0.028959</td><td align="right">+8.81</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">-2.54</td><td align="right">19</td><td align="right">1</td><td align="right">7</td><td align="right">0.028959</td><td align="right">-1.70</td></tr>
</tbody>
</table>

## Without rbg

Instances: 24

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. delta [pp]</th><th>Wins</th><th>Ties</th><th>Losses</th><th>Sign-test p-value</th><th>Success delta [pp]</th></tr>
</thead>
<tbody>
<tr><td>Strict MSA</td><td align="right">-0.19</td><td align="right">13</td><td align="right">3</td><td align="right">8</td><td align="right">0.383310</td><td align="right">+1.00</td></tr>
<tr><td>Rooted MSA</td><td align="right">+0.02</td><td align="right">13</td><td align="right">3</td><td align="right">8</td><td align="right">0.383310</td><td align="right">+1.42</td></tr>
<tr><td>Cycle cover</td><td align="right">-0.08</td><td align="right">14</td><td align="right">3</td><td align="right">7</td><td align="right">0.189247</td><td align="right">-2.42</td></tr>
<tr><td>Cycle-cover patching</td><td align="right">-0.66</td><td align="right">16</td><td align="right">1</td><td align="right">7</td><td align="right">0.093140</td><td align="right">-2.58</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">-0.78</td><td align="right">16</td><td align="right">1</td><td align="right">7</td><td align="right">0.093140</td><td align="right">-2.33</td></tr>
</tbody>
</table>

