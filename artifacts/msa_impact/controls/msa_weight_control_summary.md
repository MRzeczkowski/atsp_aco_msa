# MSA Weight Control Summary

This report checks whether MSA keeps beating control heuristics across the tested MSA heuristic-strength sweep. Each row compares MSA and a control at the same heuristic weight. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.

## Random Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Random sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.18</td><td align="right">2.58</td><td align="right">-0.41</td><td align="right">17.92</td><td align="right">15.00</td></tr>
<tr><td align="right">0.20</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">1.99</td><td align="right">2.81</td><td align="right">-0.82</td><td align="right">20.42</td><td align="right">12.08</td></tr>
<tr><td align="right">0.30</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.80</td><td align="right">3.20</td><td align="right">-1.40</td><td align="right">18.75</td><td align="right">12.08</td></tr>
<tr><td align="right">0.40</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.81</td><td align="right">3.60</td><td align="right">-1.79</td><td align="right">18.33</td><td align="right">10.97</td></tr>
<tr><td align="right">0.50</td><td align="right">8</td><td align="right">7</td><td align="right">1</td><td align="right">0</td><td align="right">0.070312</td><td align="right">2.03</td><td align="right">4.00</td><td align="right">-1.98</td><td align="right">18.33</td><td align="right">10.42</td></tr>
<tr><td align="right">0.60</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.06</td><td align="right">4.32</td><td align="right">-2.26</td><td align="right">18.75</td><td align="right">11.11</td></tr>
<tr><td align="right">0.70</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.93</td><td align="right">4.58</td><td align="right">-2.64</td><td align="right">20.00</td><td align="right">10.28</td></tr>
<tr><td align="right">0.80</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.08</td><td align="right">4.80</td><td align="right">-2.72</td><td align="right">18.33</td><td align="right">10.00</td></tr>
<tr><td align="right">0.90</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.00</td><td align="right">5.11</td><td align="right">-3.11</td><td align="right">18.75</td><td align="right">10.56</td></tr>
<tr><td align="right">1.00</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.08</td><td align="right">5.30</td><td align="right">-3.22</td><td align="right">17.92</td><td align="right">9.86</td></tr>
</tbody>
</table>

## Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.18</td><td align="right">2.36</td><td align="right">-0.18</td><td align="right">17.92</td><td align="right">14.58</td></tr>
<tr><td align="right">0.20</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.99</td><td align="right">2.52</td><td align="right">-0.53</td><td align="right">20.42</td><td align="right">15.00</td></tr>
<tr><td align="right">0.30</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.80</td><td align="right">2.45</td><td align="right">-0.65</td><td align="right">18.75</td><td align="right">14.58</td></tr>
<tr><td align="right">0.40</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.81</td><td align="right">2.54</td><td align="right">-0.73</td><td align="right">18.33</td><td align="right">15.42</td></tr>
<tr><td align="right">0.50</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.03</td><td align="right">2.51</td><td align="right">-0.48</td><td align="right">18.33</td><td align="right">17.92</td></tr>
<tr><td align="right">0.60</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.06</td><td align="right">2.76</td><td align="right">-0.69</td><td align="right">18.75</td><td align="right">14.17</td></tr>
<tr><td align="right">0.70</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.93</td><td align="right">2.61</td><td align="right">-0.67</td><td align="right">20.00</td><td align="right">17.50</td></tr>
<tr><td align="right">0.80</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.08</td><td align="right">2.60</td><td align="right">-0.52</td><td align="right">18.33</td><td align="right">17.50</td></tr>
<tr><td align="right">0.90</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.00</td><td align="right">2.65</td><td align="right">-0.65</td><td align="right">18.75</td><td align="right">15.00</td></tr>
<tr><td align="right">1.00</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.08</td><td align="right">2.68</td><td align="right">-0.60</td><td align="right">17.92</td><td align="right">17.50</td></tr>
</tbody>
</table>

## Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.18</td><td align="right">2.45</td><td align="right">-0.27</td><td align="right">17.92</td><td align="right">15.56</td></tr>
<tr><td align="right">0.20</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.99</td><td align="right">2.79</td><td align="right">-0.81</td><td align="right">20.42</td><td align="right">13.33</td></tr>
<tr><td align="right">0.30</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.80</td><td align="right">3.19</td><td align="right">-1.39</td><td align="right">18.75</td><td align="right">11.81</td></tr>
<tr><td align="right">0.40</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.81</td><td align="right">3.59</td><td align="right">-1.78</td><td align="right">18.33</td><td align="right">11.25</td></tr>
<tr><td align="right">0.50</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.03</td><td align="right">3.96</td><td align="right">-1.93</td><td align="right">18.33</td><td align="right">9.44</td></tr>
<tr><td align="right">0.60</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.06</td><td align="right">4.26</td><td align="right">-2.20</td><td align="right">18.75</td><td align="right">9.86</td></tr>
<tr><td align="right">0.70</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.93</td><td align="right">4.59</td><td align="right">-2.65</td><td align="right">20.00</td><td align="right">9.03</td></tr>
<tr><td align="right">0.80</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.08</td><td align="right">4.82</td><td align="right">-2.74</td><td align="right">18.33</td><td align="right">8.75</td></tr>
<tr><td align="right">0.90</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.00</td><td align="right">5.11</td><td align="right">-3.11</td><td align="right">18.75</td><td align="right">8.33</td></tr>
<tr><td align="right">1.00</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.08</td><td align="right">5.36</td><td align="right">-3.27</td><td align="right">17.92</td><td align="right">8.89</td></tr>
</tbody>
</table>

