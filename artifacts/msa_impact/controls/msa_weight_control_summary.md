# MSA Weight Control Summary

This report checks whether MSA keeps beating the cheaper control heuristics at the best MSA weight selected for each instance. Rows are grouped by weight because different instances may select different best weights. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.

## Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.20</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.83</td><td align="right">3.74</td><td align="right">-0.91</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.30</td><td align="right">3</td><td align="right">3</td><td align="right">0</td><td align="right">0</td><td align="right">0.250000</td><td align="right">1.94</td><td align="right">2.30</td><td align="right">-0.36</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.50</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.70</td><td align="right">2.14</td><td align="right">-1.44</td><td align="right">40.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.80</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.74</td><td align="right">1.55</td><td align="right">-0.81</td><td align="right">0.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

## Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.20</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.83</td><td align="right">6.69</td><td align="right">-3.86</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.30</td><td align="right">3</td><td align="right">3</td><td align="right">0</td><td align="right">0</td><td align="right">0.250000</td><td align="right">1.94</td><td align="right">3.00</td><td align="right">-1.06</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.50</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.70</td><td align="right">4.47</td><td align="right">-3.77</td><td align="right">40.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.80</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.74</td><td align="right">3.51</td><td align="right">-2.77</td><td align="right">0.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

