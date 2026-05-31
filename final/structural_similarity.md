# Structural Similarity To Found Optimal Tours

This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.

Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.

## Findings

- **Precision vs found-optimal tours: MSA heuristic 65.44%, cycle cover 77.15%, cycle-cover MSA patching 77.40%.**
- **Recall vs found-optimal tours: MSA heuristic 3.18%, cycle cover 6.02%, cycle-cover MSA patching 6.04%.**
- **Best-or-tied precision counts: MSA heuristic 6/23, cycle cover 6/23, cycle-cover MSA patching 13/23.**
- **Best-or-tied recall counts: MSA heuristic 0/23, cycle cover 8/23, cycle-cover MSA patching 15/23.**

<table>
<thead>
<tr><th rowspan="2">Instance</th><th colspan="2">MSA heuristic</th><th colspan="2">Cycle cover</th><th colspan="2">Cycle-cover MSA patching</th></tr>
<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>atex5</td><td align="right"><strong>100.00</strong></td><td align="right">9.07</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>20.40</strong></td><td align="right">77.78</td><td align="right">15.86</td></tr>
<tr><td>code198</td><td align="right">37.97</td><td align="right">0.26</td><td align="right">97.98</td><td align="right">0.71</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>0.73</strong></td></tr>
<tr><td>crane100_0</td><td align="right">53.57</td><td align="right">30.00</td><td align="right">53.00</td><td align="right">53.00</td><td align="right"><strong>60.00</strong></td><td align="right"><strong>60.00</strong></td></tr>
<tr><td>crane100_1</td><td align="right">57.89</td><td align="right">33.00</td><td align="right">60.00</td><td align="right">60.00</td><td align="right"><strong>65.00</strong></td><td align="right"><strong>65.00</strong></td></tr>
<tr><td>crane100_2</td><td align="right"><strong>63.33</strong></td><td align="right">38.00</td><td align="right">58.00</td><td align="right"><strong>58.00</strong></td><td align="right">55.00</td><td align="right">55.00</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>61.11</strong></td><td align="right">33.33</td><td align="right">40.91</td><td align="right">40.91</td><td align="right">46.97</td><td align="right"><strong>46.97</strong></td></tr>
<tr><td>crane66_1</td><td align="right"><strong>54.55</strong></td><td align="right">27.27</td><td align="right"><strong>54.55</strong></td><td align="right"><strong>54.55</strong></td><td align="right">37.88</td><td align="right">37.88</td></tr>
<tr><td>crane66_2</td><td align="right"><strong>73.68</strong></td><td align="right">42.42</td><td align="right">48.48</td><td align="right">48.48</td><td align="right">53.03</td><td align="right"><strong>53.03</strong></td></tr>
<tr><td>ft53</td><td align="right">53.57</td><td align="right">28.30</td><td align="right">64.15</td><td align="right">64.15</td><td align="right"><strong>69.81</strong></td><td align="right"><strong>69.81</strong></td></tr>
<tr><td>ft70</td><td align="right">51.92</td><td align="right">38.57</td><td align="right"><strong>65.71</strong></td><td align="right"><strong>65.71</strong></td><td align="right">57.14</td><td align="right">57.14</td></tr>
<tr><td>ftv100</td><td align="right">74.14</td><td align="right">31.85</td><td align="right"><strong>86.14</strong></td><td align="right"><strong>64.44</strong></td><td align="right">80.20</td><td align="right">60.00</td></tr>
<tr><td>ftv110</td><td align="right">68.75</td><td align="right">30.99</td><td align="right"><strong>79.28</strong></td><td align="right"><strong>61.97</strong></td><td align="right">78.38</td><td align="right">61.27</td></tr>
<tr><td>ftv120</td><td align="right">68.49</td><td align="right">31.65</td><td align="right">76.86</td><td align="right">58.86</td><td align="right"><strong>80.17</strong></td><td align="right"><strong>61.39</strong></td></tr>
<tr><td>ftv130</td><td align="right">75.00</td><td align="right">34.68</td><td align="right">80.92</td><td align="right">61.27</td><td align="right"><strong>82.44</strong></td><td align="right"><strong>62.43</strong></td></tr>
<tr><td>ftv140</td><td align="right">75.29</td><td align="right">34.97</td><td align="right">80.85</td><td align="right">62.30</td><td align="right"><strong>81.56</strong></td><td align="right"><strong>62.84</strong></td></tr>
<tr><td>ftv150</td><td align="right">76.09</td><td align="right">34.83</td><td align="right">84.11</td><td align="right">63.18</td><td align="right"><strong>84.77</strong></td><td align="right"><strong>63.68</strong></td></tr>
<tr><td>ftv160</td><td align="right">76.60</td><td align="right">33.03</td><td align="right">88.82</td><td align="right">65.60</td><td align="right"><strong>91.30</strong></td><td align="right"><strong>67.43</strong></td></tr>
<tr><td>ftv170</td><td align="right">76.24</td><td align="right">33.19</td><td align="right">82.46</td><td align="right">60.78</td><td align="right"><strong>83.04</strong></td><td align="right"><strong>61.21</strong></td></tr>
<tr><td>ftv55</td><td align="right"><strong>80.00</strong></td><td align="right">30.00</td><td align="right">73.21</td><td align="right"><strong>51.25</strong></td><td align="right">71.43</td><td align="right">50.00</td></tr>
<tr><td>ftv64</td><td align="right">76.32</td><td align="right">31.52</td><td align="right">73.85</td><td align="right">52.17</td><td align="right"><strong>86.15</strong></td><td align="right"><strong>60.87</strong></td></tr>
<tr><td>ftv70</td><td align="right">71.43</td><td align="right">32.26</td><td align="right">70.42</td><td align="right">53.76</td><td align="right"><strong>77.46</strong></td><td align="right"><strong>59.14</strong></td></tr>
<tr><td>ftv90</td><td align="right">73.08</td><td align="right">33.63</td><td align="right"><strong>82.42</strong></td><td align="right"><strong>66.37</strong></td><td align="right">76.92</td><td align="right">61.95</td></tr>
<tr><td>td100_1</td><td align="right">57.32</td><td align="right">16.79</td><td align="right">97.03</td><td align="right">35.00</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>36.07</strong></td></tr>
<tr><td><strong>Total</strong></td><td align="right">65.44</td><td align="right">3.18</td><td align="right">77.15</td><td align="right">6.02</td><td align="right"><strong>77.40</strong></td><td align="right"><strong>6.04</strong></td></tr>
</tbody>
</table>
