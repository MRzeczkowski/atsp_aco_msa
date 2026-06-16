# MSA Weight Control Summary

This report checks whether each MSA-impact variant keeps beating its matched control heuristics at the best variant weight selected for each instance. Rows are grouped by weight because different instances may select different best weights. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.

## Strict MSA vs Strict Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>Strict MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>Strict MSA avg dev. [%]</th><th>Strict distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>Strict MSA success [%]</th><th>Strict distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.03</td><td align="right">1.13</td><td align="right">-0.10</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.40</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.79</td><td align="right">4.42</td><td align="right">-1.63</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.50</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.36</td><td align="right">0.72</td><td align="right">-0.36</td><td align="right">50.00</td><td align="right">20.00</td></tr>
<tr><td align="right">0.80</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">1.31</td><td align="right">2.65</td><td align="right">-1.33</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.12</td><td align="right">3.47</td><td align="right">-1.35</td><td align="right">20.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

## Strict MSA vs Strict Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>Strict MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>Strict MSA avg dev. [%]</th><th>Strict shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>Strict MSA success [%]</th><th>Strict shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.03</td><td align="right">1.49</td><td align="right">-0.46</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.40</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.79</td><td align="right">7.43</td><td align="right">-4.64</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.50</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.36</td><td align="right">2.21</td><td align="right">-1.85</td><td align="right">50.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.80</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">1.31</td><td align="right">5.43</td><td align="right">-4.11</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.12</td><td align="right">4.74</td><td align="right">-2.62</td><td align="right">20.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

## Rooted MSA vs Rooted Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>Rooted MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>Rooted MSA avg dev. [%]</th><th>Rooted distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>Rooted MSA success [%]</th><th>Rooted distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.97</td><td align="right">3.85</td><td align="right">-0.88</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.05</td><td align="right">1.20</td><td align="right">-0.15</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.30</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.45</td><td align="right">1.24</td><td align="right">-0.79</td><td align="right">40.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.40</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.76</td><td align="right">3.02</td><td align="right">-2.26</td><td align="right">10.00</td><td align="right">10.00</td></tr>
<tr><td align="right">0.80</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.70</td><td align="right">1.21</td><td align="right">-0.51</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.67</td><td align="right">4.10</td><td align="right">-1.43</td><td align="right">0.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

## Rooted MSA vs Rooted Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>Rooted MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>Rooted MSA avg dev. [%]</th><th>Rooted shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>Rooted MSA success [%]</th><th>Rooted shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">2</td><td align="right">2</td><td align="right">0</td><td align="right">0</td><td align="right">0.500000</td><td align="right">2.97</td><td align="right">4.54</td><td align="right">-1.56</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.20</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">1.05</td><td align="right">1.40</td><td align="right">-0.35</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.30</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.45</td><td align="right">1.23</td><td align="right">-0.78</td><td align="right">40.00</td><td align="right">23.33</td></tr>
<tr><td align="right">0.40</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.76</td><td align="right">3.61</td><td align="right">-2.85</td><td align="right">10.00</td><td align="right">0.00</td></tr>
<tr><td align="right">0.80</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">0.70</td><td align="right">2.83</td><td align="right">-2.13</td><td align="right">0.00</td><td align="right">0.00</td></tr>
<tr><td align="right">1.00</td><td align="right">1</td><td align="right">1</td><td align="right">0</td><td align="right">0</td><td align="right">1.000000</td><td align="right">2.67</td><td align="right">11.46</td><td align="right">-8.79</td><td align="right">0.00</td><td align="right">0.00</td></tr>
</tbody>
</table>

