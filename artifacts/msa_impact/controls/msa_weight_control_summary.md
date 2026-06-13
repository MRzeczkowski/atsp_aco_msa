# MSA Weight Control Summary

This report checks whether MSA keeps beating the cheaper control heuristics at the best MSA weight selected for each instance. Rows are grouped by weight because different instances may select different best weights. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.

## Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.61</td><td align="right">3.59</td><td align="right">-0.98</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.06</td><td align="right">1.13</td><td align="right">-0.07</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.70</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.36</td><td align="right">1.12</td><td align="right">-0.76</td><td align="right">60.00</td><td align="right">10.00</td></tr>
<tr><td align="right">0.90</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">1.85</td><td align="right">3.20</td><td align="right">-1.34</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.10</td><td align="right">4.78</td><td align="right">-2.67</td><td align="right">10.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

## Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.61</td><td align="right">4.41</td><td align="right">-1.80</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.06</td><td align="right">1.47</td><td align="right">-0.41</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.70</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.36</td><td align="right">2.85</td><td align="right">-2.49</td><td align="right">60.00</td><td align="right">3.33</td></tr>
<tr><td align="right">0.90</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">1.85</td><td align="right">6.12</td><td align="right">-4.26</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.10</td><td align="right">6.90</td><td align="right">-4.80</td><td align="right">10.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

