# MSA Weight Control Summary

This report checks whether MSA keeps beating the cheaper control heuristics at the best MSA weight selected for each instance. Rows are grouped by weight because different instances may select different best weights. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.

## Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.03</td><td align="right">1.13</td><td align="right">-0.10</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.40</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.79</td><td align="right">4.42</td><td align="right">-1.63</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.50</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.36</td><td align="right">0.72</td><td align="right">-0.36</td><td align="right">50.00</td><td align="right">20.00</td></tr>
<tr><td align="right">0.80</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">1.31</td><td align="right">2.65</td><td align="right">-1.33</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.12</td><td align="right">3.47</td><td align="right">-1.35</td><td align="right">20.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

## Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.03</td><td align="right">1.49</td><td align="right">-0.46</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.40</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.79</td><td align="right">7.43</td><td align="right">-4.64</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.50</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.36</td><td align="right">2.21</td><td align="right">-1.85</td><td align="right">50.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.80</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">1.31</td><td align="right">5.43</td><td align="right">-4.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.12</td><td align="right">4.74</td><td align="right">-2.62</td><td align="right">20.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

