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
<tr><td>Strict MSA</td><td align="right">-0.02</td><td align="right">14</td><td align="right">8</td><td align="right">5</td><td align="right">0.063568</td><td align="right">+2.67</td></tr>
<tr><td>Rooted MSA</td><td align="right">-0.02</td><td align="right">16</td><td align="right">5</td><td align="right">6</td><td align="right">0.052479</td><td align="right">+4.44</td></tr>
<tr><td>Cycle cover</td><td align="right">-0.28</td><td align="right">13</td><td align="right">3</td><td align="right">11</td><td align="right">0.838820</td><td align="right">-5.70</td></tr>
<tr><td>Cycle-cover patching</td><td align="right">-0.26</td><td align="right">14</td><td align="right">2</td><td align="right">11</td><td align="right">0.690038</td><td align="right">+8.07</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">-0.26</td><td align="right">11</td><td align="right">2</td><td align="right">14</td><td align="right">0.690038</td><td align="right">-4.15</td></tr>
</tbody>
</table>

## Without rbg

Instances: 24

<table>
<thead>
<tr><th>Heuristic</th><th>Avg best dev. delta [pp]</th><th>Wins</th><th>Ties</th><th>Losses</th><th>Sign-test p-value</th><th>Success delta [pp]</th></tr>
</thead>
<tbody>
<tr><td>Strict MSA</td><td align="right">-0.02</td><td align="right">12</td><td align="right">8</td><td align="right">4</td><td align="right">0.076813</td><td align="right">+3.00</td></tr>
<tr><td>Rooted MSA</td><td align="right">-0.03</td><td align="right">14</td><td align="right">5</td><td align="right">5</td><td align="right">0.063568</td><td align="right">+5.00</td></tr>
<tr><td>Cycle cover</td><td align="right">-0.07</td><td align="right">10</td><td align="right">3</td><td align="right">11</td><td align="right">1.000000</td><td align="right">-7.17</td></tr>
<tr><td>Cycle-cover patching</td><td align="right">-0.03</td><td align="right">11</td><td align="right">2</td><td align="right">11</td><td align="right">1.000000</td><td align="right">-3.42</td></tr>
<tr><td>Cycle-cover MSA patching</td><td align="right">-0.04</td><td align="right">8</td><td align="right">2</td><td align="right">14</td><td align="right">0.286279</td><td align="right">-8.42</td></tr>
</tbody>
</table>

