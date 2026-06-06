# Structural Similarity To Found Optimal Tours

This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.

Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.

## Findings

- **Precision vs found-optimal tours: MSA heuristic 61.39%, cycle cover 79.88%, cycle-cover MSA patching 77.87%.**
- **Recall vs found-optimal tours: MSA heuristic 1.55%, cycle cover 3.15%, cycle-cover MSA patching 3.07%.**
- **Best-or-tied precision counts: MSA heuristic 5/13, cycle cover 3/13, cycle-cover MSA patching 7/13.**
- **Best-or-tied recall counts: MSA heuristic 0/13, cycle cover 6/13, cycle-cover MSA patching 7/13.**

<table>
<thead>
<tr><th rowspan="2">Instance</th><th colspan="2">MSA heuristic</th><th colspan="2">Cycle cover</th><th colspan="2">Cycle-cover MSA patching</th></tr>
<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>100.00</strong></td><td align="right">8.99</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>20.22</strong></td><td align="right">81.94</td><td align="right">16.57</td></tr>
<tr><td>code198</td><td align="right">37.97</td><td align="right">0.26</td><td align="right">97.98</td><td align="right">0.71</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>0.73</strong></td></tr>
<tr><td>crane100_1</td><td align="right">57.89</td><td align="right">33.00</td><td align="right">60.00</td><td align="right">60.00</td><td align="right"><strong>65.00</strong></td><td align="right"><strong>65.00</strong></td></tr>
<tr><td>crane66_1</td><td align="right"><strong>54.55</strong></td><td align="right">27.27</td><td align="right"><strong>54.55</strong></td><td align="right"><strong>54.55</strong></td><td align="right">37.88</td><td align="right">37.88</td></tr>
<tr><td>ft53</td><td align="right">53.57</td><td align="right">28.30</td><td align="right">64.15</td><td align="right">64.15</td><td align="right"><strong>69.81</strong></td><td align="right"><strong>69.81</strong></td></tr>
<tr><td>ftv120</td><td align="right">68.49</td><td align="right">31.65</td><td align="right">76.86</td><td align="right">58.86</td><td align="right"><strong>80.17</strong></td><td align="right"><strong>61.39</strong></td></tr>
<tr><td>ftv150</td><td align="right">76.09</td><td align="right">34.83</td><td align="right">84.11</td><td align="right">63.18</td><td align="right"><strong>84.77</strong></td><td align="right"><strong>63.68</strong></td></tr>
<tr><td>ftv33</td><td align="right"><strong>70.00</strong></td><td align="right">37.84</td><td align="right">52.94</td><td align="right"><strong>48.65</strong></td><td align="right">38.24</td><td align="right">35.14</td></tr>
<tr><td>ftv64</td><td align="right">76.32</td><td align="right">31.52</td><td align="right">73.85</td><td align="right">52.17</td><td align="right"><strong>86.15</strong></td><td align="right"><strong>60.87</strong></td></tr>
<tr><td>ftv90</td><td align="right">73.08</td><td align="right">33.63</td><td align="right"><strong>82.42</strong></td><td align="right"><strong>66.37</strong></td><td align="right">76.92</td><td align="right">61.95</td></tr>
<tr><td>p43</td><td align="right"><strong>100.00</strong></td><td align="right">12.87</td><td align="right">88.37</td><td align="right"><strong>22.22</strong></td><td align="right">60.47</td><td align="right">15.20</td></tr>
<tr><td>ry48p</td><td align="right"><strong>64.71</strong></td><td align="right">22.92</td><td align="right">41.67</td><td align="right"><strong>41.67</strong></td><td align="right">31.25</td><td align="right">31.25</td></tr>
<tr><td>td100_1</td><td align="right">57.32</td><td align="right">16.79</td><td align="right">97.03</td><td align="right">35.00</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>36.07</strong></td></tr>
<tr><td><strong>Total</strong></td><td align="right">61.39</td><td align="right">1.55</td><td align="right"><strong>79.88</strong></td><td align="right"><strong>3.15</strong></td><td align="right">77.87</td><td align="right">3.07</td></tr>
</tbody>
</table>
