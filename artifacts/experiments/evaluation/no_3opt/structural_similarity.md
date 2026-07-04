# Structural Similarity To Found Optimal Tours

This table compares the current MSA heuristic, minimum cycle-cover, GKS patching, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.

Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.

## Findings

- **Precision vs found-optimal tours: MSA heuristic 65.07%, cycle cover 70.61%, GKS patching 71.54%, cycle-cover MSA patching 70.34%.**
- **Recall vs found-optimal tours: MSA heuristic 2.90%, cycle cover 4.42%, GKS patching 4.48%, cycle-cover MSA patching 4.41%.**
- **Best-or-tied precision counts: MSA heuristic 10/24, cycle cover 5/24, GKS patching 10/24, cycle-cover MSA patching 5/24.**
- **Best-or-tied recall counts: MSA heuristic 0/24, cycle cover 8/24, GKS patching 13/24, cycle-cover MSA patching 10/24.**

<table>
<thead>
<tr><th rowspan="2">Instance</th><th colspan="2">MSA heuristic</th><th colspan="2">Cycle cover</th><th colspan="2">GKS patching</th><th colspan="2">Cycle-cover MSA patching</th></tr>
<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>
</thead>
<tbody>
<tr><td>atex1</td><td align="right"><strong>57.14</strong></td><td align="right">17.39</td><td align="right">50.00</td><td align="right"><strong>34.78</strong></td><td align="right">50.00</td><td align="right"><strong>34.78</strong></td><td align="right">50.00</td><td align="right"><strong>34.78</strong></td></tr>
<tr><td>atex3</td><td align="right"><strong>100.00</strong></td><td align="right">12.50</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>25.00</strong></td><td align="right"><strong>100.00</strong></td><td align="right"><strong>25.00</strong></td><td align="right">71.88</td><td align="right">17.97</td></tr>
<tr><td>atex4</td><td align="right"><strong>100.00</strong></td><td align="right">10.65</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>18.25</strong></td><td align="right">75.00</td><td align="right">13.69</td><td align="right">70.83</td><td align="right">12.93</td></tr>
<tr><td>code198</td><td align="right">72.73</td><td align="right">0.49</td><td align="right">97.98</td><td align="right">0.70</td><td align="right"><strong>100.00</strong></td><td align="right"><strong>0.72</strong></td><td align="right"><strong>100.00</strong></td><td align="right"><strong>0.72</strong></td></tr>
<tr><td>crane100_0</td><td align="right">54.55</td><td align="right">30.00</td><td align="right">53.00</td><td align="right">53.00</td><td align="right"><strong>64.00</strong></td><td align="right"><strong>64.00</strong></td><td align="right"><strong>64.00</strong></td><td align="right"><strong>64.00</strong></td></tr>
<tr><td>crane100_2</td><td align="right"><strong>63.33</strong></td><td align="right">38.00</td><td align="right">58.00</td><td align="right"><strong>58.00</strong></td><td align="right">56.00</td><td align="right">56.00</td><td align="right">56.00</td><td align="right">56.00</td></tr>
<tr><td>crane66_0</td><td align="right"><strong>61.11</strong></td><td align="right">33.33</td><td align="right">40.91</td><td align="right"><strong>40.91</strong></td><td align="right">34.85</td><td align="right">34.85</td><td align="right">40.91</td><td align="right"><strong>40.91</strong></td></tr>
<tr><td>crane66_2</td><td align="right"><strong>76.32</strong></td><td align="right">43.94</td><td align="right">48.48</td><td align="right">48.48</td><td align="right">46.97</td><td align="right">46.97</td><td align="right">53.03</td><td align="right"><strong>53.03</strong></td></tr>
<tr><td>ft70</td><td align="right">52.83</td><td align="right">40.00</td><td align="right"><strong>65.71</strong></td><td align="right"><strong>65.71</strong></td><td align="right"><strong>65.71</strong></td><td align="right"><strong>65.71</strong></td><td align="right">61.43</td><td align="right">61.43</td></tr>
<tr><td>ftv100</td><td align="right">74.14</td><td align="right">31.85</td><td align="right"><strong>86.14</strong></td><td align="right"><strong>64.44</strong></td><td align="right">84.16</td><td align="right">62.96</td><td align="right">80.20</td><td align="right">60.00</td></tr>
<tr><td>ftv110</td><td align="right">68.75</td><td align="right">30.99</td><td align="right"><strong>79.28</strong></td><td align="right"><strong>61.97</strong></td><td align="right">78.38</td><td align="right">61.27</td><td align="right">75.68</td><td align="right">59.15</td></tr>
<tr><td>ftv130</td><td align="right">75.00</td><td align="right">34.68</td><td align="right">80.92</td><td align="right">61.27</td><td align="right"><strong>82.44</strong></td><td align="right"><strong>62.43</strong></td><td align="right">80.15</td><td align="right">60.69</td></tr>
<tr><td>ftv140</td><td align="right">75.29</td><td align="right">34.97</td><td align="right">80.85</td><td align="right">62.30</td><td align="right"><strong>83.69</strong></td><td align="right"><strong>64.48</strong></td><td align="right">79.43</td><td align="right">61.20</td></tr>
<tr><td>ftv160</td><td align="right">76.60</td><td align="right">33.03</td><td align="right">88.82</td><td align="right">65.60</td><td align="right"><strong>94.41</strong></td><td align="right"><strong>69.72</strong></td><td align="right">91.30</td><td align="right">67.43</td></tr>
<tr><td>ftv170</td><td align="right">76.24</td><td align="right">33.19</td><td align="right">82.46</td><td align="right">60.78</td><td align="right"><strong>88.89</strong></td><td align="right"><strong>65.52</strong></td><td align="right">81.87</td><td align="right">60.34</td></tr>
<tr><td>ftv35</td><td align="right"><strong>61.90</strong></td><td align="right">36.11</td><td align="right">44.44</td><td align="right">44.44</td><td align="right">36.11</td><td align="right">36.11</td><td align="right">47.22</td><td align="right"><strong>47.22</strong></td></tr>
<tr><td>ftv38</td><td align="right"><strong>60.87</strong></td><td align="right">33.33</td><td align="right">46.15</td><td align="right">42.86</td><td align="right">38.46</td><td align="right">35.71</td><td align="right">48.72</td><td align="right"><strong>45.24</strong></td></tr>
<tr><td>ftv44</td><td align="right"><strong>69.23</strong></td><td align="right">31.03</td><td align="right">51.11</td><td align="right">39.66</td><td align="right">53.33</td><td align="right"><strong>41.38</strong></td><td align="right">48.89</td><td align="right">37.93</td></tr>
<tr><td>ftv47</td><td align="right">60.00</td><td align="right">25.00</td><td align="right">62.50</td><td align="right">50.00</td><td align="right">72.92</td><td align="right">58.33</td><td align="right"><strong>75.00</strong></td><td align="right"><strong>60.00</strong></td></tr>
<tr><td>ftv55</td><td align="right"><strong>80.65</strong></td><td align="right">31.25</td><td align="right">73.21</td><td align="right">51.25</td><td align="right">76.79</td><td align="right"><strong>53.75</strong></td><td align="right">71.43</td><td align="right">50.00</td></tr>
<tr><td>ftv70</td><td align="right">71.43</td><td align="right">32.26</td><td align="right">70.42</td><td align="right">53.76</td><td align="right">73.24</td><td align="right">55.91</td><td align="right"><strong>77.46</strong></td><td align="right"><strong>59.14</strong></td></tr>
<tr><td>rbg358</td><td align="right">56.39</td><td align="right">7.35</td><td align="right">63.41</td><td align="right">9.70</td><td align="right">63.41</td><td align="right">9.70</td><td align="right"><strong>63.97</strong></td><td align="right"><strong>9.79</strong></td></tr>
<tr><td>rbg403</td><td align="right">57.57</td><td align="right">2.54</td><td align="right">65.76</td><td align="right">3.47</td><td align="right"><strong>66.25</strong></td><td align="right"><strong>3.49</strong></td><td align="right">65.76</td><td align="right">3.47</td></tr>
<tr><td>rbg443</td><td align="right">59.95</td><td align="right">2.69</td><td align="right">62.98</td><td align="right">3.37</td><td align="right"><strong>63.66</strong></td><td align="right"><strong>3.40</strong></td><td align="right">62.75</td><td align="right">3.36</td></tr>
<tr><td><strong>Total</strong></td><td align="right">65.07</td><td align="right">2.90</td><td align="right">70.61</td><td align="right">4.42</td><td align="right"><strong>71.54</strong></td><td align="right"><strong>4.48</strong></td><td align="right">70.34</td><td align="right">4.41</td></tr>
</tbody>
</table>
