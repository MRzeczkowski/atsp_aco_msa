# Structural Similarity To Found Optimal Tours

This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.

Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.

## Findings

- **Precision vs found-optimal tours: MSA heuristic 75.13%, cycle cover 73.06%, cycle-cover MSA patching 57.78%.**
- **Recall vs found-optimal tours: MSA heuristic 20.98%, cycle cover 38.06%, cycle-cover MSA patching 30.10%.**
- **Best-or-tied precision counts: MSA heuristic 4/5, cycle cover 3/5, cycle-cover MSA patching 0/5.**
- **Best-or-tied recall counts: MSA heuristic 0/5, cycle cover 5/5, cycle-cover MSA patching 0/5.**

<table>
<thead>
<tr><th rowspan="2">Instance</th><th colspan="2">MSA heuristic</th><th colspan="2">Cycle cover</th><th colspan="2">Cycle-cover MSA patching</th></tr>
<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>100.00</strong></td><td align="right">9.04</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>20.34</strong></td><td align="right">61.11</td><td align="right">12.43</td></tr>
<tr><td>crane66_1</td><td align="right"><strong>54.55</strong></td><td align="right">27.27</td><td align="right"><strong>54.55</strong></td><td align="right"><strong>54.55</strong></td><td align="right">34.85</td><td align="right">34.85</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>73.68</strong></td><td align="right">42.42</td><td align="right">48.48</td><td align="right"><strong>48.48</strong></td><td align="right">40.91</td><td align="right">40.91</td></tr>
<tr><td>ftv64</td><td align="right"><strong>76.32</strong></td><td align="right">31.52</td><td align="right">73.85</td><td align="right"><strong>52.17</strong></td><td align="right">67.69</td><td align="right">47.83</td></tr>
<tr><td>ftv90</td><td align="right">73.08</td><td align="right">33.63</td><td align="right"><strong>82.42</strong></td><td align="right"><strong>66.37</strong></td><td align="right">76.92</td><td align="right">61.95</td></tr>
<tr><td><strong>Total</strong></td><td align="right"><strong>75.13</strong></td><td align="right">20.98</td><td align="right">73.06</td><td align="right"><strong>38.06</strong></td><td align="right">57.78</td><td align="right">30.10</td></tr>
</tbody>
</table>
