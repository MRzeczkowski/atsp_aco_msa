# Structural Similarity To Found Optimal Tours

This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.

Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.

## Findings

- **Precision vs found-optimal tours: MSA heuristic 67.13%, cycle cover 78.20%, cycle-cover MSA patching 76.73%.**
- **Recall vs found-optimal tours: MSA heuristic 19.34%, cycle cover 37.56%, cycle-cover MSA patching 36.86%.**
- **Best-or-tied precision counts: MSA heuristic 2/6, cycle cover 3/6, cycle-cover MSA patching 3/6.**
- **Best-or-tied recall counts: MSA heuristic 0/6, cycle cover 4/6, cycle-cover MSA patching 3/6.**

<table>
<thead>
<tr><th rowspan="2">Instance</th><th colspan="2">MSA heuristic</th><th colspan="2">Cycle cover</th><th colspan="2">Cycle-cover MSA patching</th></tr>
<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>100.00</strong></td><td align="right">8.89</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>20.00</strong></td><td align="right">86.11</td><td align="right">17.22</td></tr>
<tr><td>crane100_1</td><td align="right">54.69</td><td align="right">35.00</td><td align="right"><strong>60.00</strong></td><td align="right"><strong>60.00</strong></td><td align="right"><strong>60.00</strong></td><td align="right"><strong>60.00</strong></td></tr>
<tr><td>ftv64</td><td align="right">76.32</td><td align="right">31.52</td><td align="right">73.85</td><td align="right">52.17</td><td align="right"><strong>86.15</strong></td><td align="right"><strong>60.87</strong></td></tr>
<tr><td>ftv90</td><td align="right">73.08</td><td align="right">33.63</td><td align="right"><strong>82.42</strong></td><td align="right"><strong>66.37</strong></td><td align="right">78.02</td><td align="right">62.83</td></tr>
<tr><td>ry48p</td><td align="right"><strong>64.71</strong></td><td align="right">22.92</td><td align="right">41.67</td><td align="right"><strong>41.67</strong></td><td align="right">33.33</td><td align="right">33.33</td></tr>
<tr><td>td100_1</td><td align="right">56.63</td><td align="right">16.79</td><td align="right">97.03</td><td align="right">35.00</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>36.07</strong></td></tr>
<tr><td><strong>Total</strong></td><td align="right">67.13</td><td align="right">19.34</td><td align="right"><strong>78.20</strong></td><td align="right"><strong>37.56</strong></td><td align="right">76.73</td><td align="right">36.86</td></tr>
</tbody>
</table>
