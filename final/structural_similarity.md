# Structural Similarity To Found Optimal Tours

This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.

Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.

## Findings

- **Precision vs found-optimal tours: MSA heuristic 65.44%, cycle cover 77.15%, cycle-cover MSA patching 74.28%.**
- **Recall vs found-optimal tours: MSA heuristic 3.18%, cycle cover 6.02%, cycle-cover MSA patching 6.21%.**
- **Best-or-tied precision counts: MSA heuristic 8/23, cycle cover 16/23, cycle-cover MSA patching 2/23.**
- **Best-or-tied recall counts: MSA heuristic 0/23, cycle cover 1/23, cycle-cover MSA patching 23/23.**

<table>
<thead>
<tr><th rowspan="2">Instance</th><th colspan="2">MSA heuristic</th><th colspan="2">Cycle cover</th><th colspan="2">Cycle-cover MSA patching</th></tr>
<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>100.00</strong></td><td align="right">9.07</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>20.40</strong></td><td align="right"><strong>100.00</strong></td><td align="right"><strong>20.40</strong></td></tr>
<tr><td>code198</td><td align="right">37.97</td><td align="right">0.26</td><td align="right"><strong>97.98</strong></td><td align="right">0.71</td><td align="right">97.50</td><td align="right"><strong>0.72</strong></td></tr>
<tr><td>crane100_0</td><td align="right"><strong>53.57</strong></td><td align="right">30.00</td><td align="right">53.00</td><td align="right">53.00</td><td align="right">51.89</td><td align="right"><strong>55.00</strong></td></tr>
<tr><td>crane100_1</td><td align="right">57.89</td><td align="right">33.00</td><td align="right"><strong>60.00</strong></td><td align="right">60.00</td><td align="right">58.72</td><td align="right"><strong>64.00</strong></td></tr>
<tr><td>crane100_2</td><td align="right"><strong>63.33</strong></td><td align="right">38.00</td><td align="right">58.00</td><td align="right">58.00</td><td align="right">55.26</td><td align="right"><strong>63.00</strong></td></tr>
<tr><td>crane66_0</td><td align="right"><strong>61.11</strong></td><td align="right">33.33</td><td align="right">40.91</td><td align="right">40.91</td><td align="right">41.43</td><td align="right"><strong>43.94</strong></td></tr>
<tr><td>crane66_1</td><td align="right"><strong>54.55</strong></td><td align="right">27.27</td><td align="right"><strong>54.55</strong></td><td align="right">54.55</td><td align="right">51.39</td><td align="right"><strong>56.06</strong></td></tr>
<tr><td>crane66_2</td><td align="right"><strong>73.68</strong></td><td align="right">42.42</td><td align="right">48.48</td><td align="right">48.48</td><td align="right">50.00</td><td align="right"><strong>56.06</strong></td></tr>
<tr><td>ft53</td><td align="right">53.57</td><td align="right">28.30</td><td align="right"><strong>64.15</strong></td><td align="right">64.15</td><td align="right">62.50</td><td align="right"><strong>66.04</strong></td></tr>
<tr><td>ft70</td><td align="right">51.92</td><td align="right">38.57</td><td align="right"><strong>65.71</strong></td><td align="right">65.71</td><td align="right">58.75</td><td align="right"><strong>67.14</strong></td></tr>
<tr><td>ftv100</td><td align="right">74.14</td><td align="right">31.85</td><td align="right"><strong>86.14</strong></td><td align="right">64.44</td><td align="right">82.41</td><td align="right"><strong>65.93</strong></td></tr>
<tr><td>ftv110</td><td align="right">68.75</td><td align="right">30.99</td><td align="right"><strong>79.28</strong></td><td align="right">61.97</td><td align="right">75.63</td><td align="right"><strong>63.38</strong></td></tr>
<tr><td>ftv120</td><td align="right">68.49</td><td align="right">31.65</td><td align="right"><strong>76.86</strong></td><td align="right">58.86</td><td align="right">73.08</td><td align="right"><strong>60.13</strong></td></tr>
<tr><td>ftv130</td><td align="right">75.00</td><td align="right">34.68</td><td align="right"><strong>80.92</strong></td><td align="right">61.27</td><td align="right">77.14</td><td align="right"><strong>62.43</strong></td></tr>
<tr><td>ftv140</td><td align="right">75.29</td><td align="right">34.97</td><td align="right"><strong>80.85</strong></td><td align="right">62.30</td><td align="right">77.48</td><td align="right"><strong>63.93</strong></td></tr>
<tr><td>ftv150</td><td align="right">76.09</td><td align="right">34.83</td><td align="right"><strong>84.11</strong></td><td align="right">63.18</td><td align="right">81.48</td><td align="right"><strong>65.67</strong></td></tr>
<tr><td>ftv160</td><td align="right">76.60</td><td align="right">33.03</td><td align="right"><strong>88.82</strong></td><td align="right">65.60</td><td align="right">86.05</td><td align="right"><strong>67.89</strong></td></tr>
<tr><td>ftv170</td><td align="right">76.24</td><td align="right">33.19</td><td align="right"><strong>82.46</strong></td><td align="right">60.78</td><td align="right">78.92</td><td align="right"><strong>62.93</strong></td></tr>
<tr><td>ftv55</td><td align="right"><strong>80.00</strong></td><td align="right">30.00</td><td align="right">73.21</td><td align="right">51.25</td><td align="right">74.14</td><td align="right"><strong>53.75</strong></td></tr>
<tr><td>ftv64</td><td align="right"><strong>76.32</strong></td><td align="right">31.52</td><td align="right">73.85</td><td align="right">52.17</td><td align="right">69.44</td><td align="right"><strong>54.35</strong></td></tr>
<tr><td>ftv70</td><td align="right">71.43</td><td align="right">32.26</td><td align="right">70.42</td><td align="right">53.76</td><td align="right"><strong>71.62</strong></td><td align="right"><strong>56.99</strong></td></tr>
<tr><td>ftv90</td><td align="right">73.08</td><td align="right">33.63</td><td align="right"><strong>82.42</strong></td><td align="right">66.37</td><td align="right">78.35</td><td align="right"><strong>67.26</strong></td></tr>
<tr><td>td100_1</td><td align="right">57.32</td><td align="right">16.79</td><td align="right"><strong>97.03</strong></td><td align="right">35.00</td><td align="right">90.00</td><td align="right"><strong>35.36</strong></td></tr>
<tr><td><strong>Total</strong></td><td align="right">65.44</td><td align="right">3.18</td><td align="right"><strong>77.15</strong></td><td align="right">6.02</td><td align="right">74.28</td><td align="right"><strong>6.21</strong></td></tr>
</tbody>
</table>
