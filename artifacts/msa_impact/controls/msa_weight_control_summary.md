# MSA Weight Control Summary

This report checks whether MSA keeps beating control heuristics across the tested MSA heuristic-strength sweep. Each row compares MSA and a control at the same heuristic weight. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.

## Random Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Random sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.16</td><td align="right">2.58</td><td align="right">-0.43</td><td align="right">17.08</td><td align="right">15.00</td></tr>
<tr><td align="right">0.20</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.17</td><td align="right">2.81</td><td align="right">-0.63</td><td align="right">17.92</td><td align="right">12.08</td></tr>
<tr><td align="right">0.30</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.89</td><td align="right">3.20</td><td align="right">-1.31</td><td align="right">17.08</td><td align="right">12.08</td></tr>
<tr><td align="right">0.40</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.80</td><td align="right">3.60</td><td align="right">-1.80</td><td align="right">17.92</td><td align="right">10.97</td></tr>
<tr><td align="right">0.50</td><td align="right">8</td><td align="right">7</td><td align="right">1</td><td align="right">0</td><td align="right">0.070312</td><td align="right">2.11</td><td align="right">4.00</td><td align="right">-1.89</td><td align="right">20.00</td><td align="right">10.42</td></tr>
<tr><td align="right">0.60</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.12</td><td align="right">4.32</td><td align="right">-2.20</td><td align="right">15.42</td><td align="right">11.11</td></tr>
<tr><td align="right">0.70</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.98</td><td align="right">4.58</td><td align="right">-2.60</td><td align="right">21.25</td><td align="right">10.28</td></tr>
<tr><td align="right">0.80</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.03</td><td align="right">4.80</td><td align="right">-2.77</td><td align="right">17.08</td><td align="right">10.00</td></tr>
<tr><td align="right">0.90</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.10</td><td align="right">5.11</td><td align="right">-3.01</td><td align="right">17.50</td><td align="right">10.56</td></tr>
<tr><td align="right">1.00</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.89</td><td align="right">5.30</td><td align="right">-3.41</td><td align="right">20.00</td><td align="right">9.86</td></tr>
</tbody>
</table>

## Distance-ranked Sparse

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Distance-ranked sparse avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked sparse success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.16</td><td align="right">2.36</td><td align="right">-0.20</td><td align="right">17.08</td><td align="right">14.58</td></tr>
<tr><td align="right">0.20</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.17</td><td align="right">2.52</td><td align="right">-0.34</td><td align="right">17.92</td><td align="right">15.00</td></tr>
<tr><td align="right">0.30</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.89</td><td align="right">2.45</td><td align="right">-0.56</td><td align="right">17.08</td><td align="right">14.58</td></tr>
<tr><td align="right">0.40</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">1.80</td><td align="right">2.54</td><td align="right">-0.74</td><td align="right">17.92</td><td align="right">15.42</td></tr>
<tr><td align="right">0.50</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.11</td><td align="right">2.51</td><td align="right">-0.39</td><td align="right">20.00</td><td align="right">17.92</td></tr>
<tr><td align="right">0.60</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.12</td><td align="right">2.76</td><td align="right">-0.64</td><td align="right">15.42</td><td align="right">14.17</td></tr>
<tr><td align="right">0.70</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.98</td><td align="right">2.61</td><td align="right">-0.63</td><td align="right">21.25</td><td align="right">17.50</td></tr>
<tr><td align="right">0.80</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.03</td><td align="right">2.60</td><td align="right">-0.57</td><td align="right">17.08</td><td align="right">17.50</td></tr>
<tr><td align="right">0.90</td><td align="right">8</td><td align="right">5</td><td align="right">2</td><td align="right">1</td><td align="right">0.453125</td><td align="right">2.10</td><td align="right">2.65</td><td align="right">-0.55</td><td align="right">17.50</td><td align="right">15.00</td></tr>
<tr><td align="right">1.00</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">1.89</td><td align="right">2.68</td><td align="right">-0.80</td><td align="right">20.00</td><td align="right">17.50</td></tr>
</tbody>
</table>

## Shuffled MSA

<table>
<thead>
<tr><th>Weight</th><th>Instances</th><th>MSA wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>MSA avg dev. [%]</th><th>Shuffled MSA avg dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Shuffled MSA success [%]</th></tr>
</thead>
<tbody>
<tr><td align="right">0.10</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.16</td><td align="right">2.50</td><td align="right">-0.34</td><td align="right">17.08</td><td align="right">14.72</td></tr>
<tr><td align="right">0.20</td><td align="right">8</td><td align="right">6</td><td align="right">1</td><td align="right">1</td><td align="right">0.125000</td><td align="right">2.17</td><td align="right">2.85</td><td align="right">-0.67</td><td align="right">17.92</td><td align="right">12.50</td></tr>
<tr><td align="right">0.30</td><td align="right">8</td><td align="right">7</td><td align="right">0</td><td align="right">1</td><td align="right">0.015625</td><td align="right">1.89</td><td align="right">3.33</td><td align="right">-1.45</td><td align="right">17.08</td><td align="right">11.94</td></tr>
<tr><td align="right">0.40</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.80</td><td align="right">3.66</td><td align="right">-1.86</td><td align="right">17.92</td><td align="right">11.25</td></tr>
<tr><td align="right">0.50</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.11</td><td align="right">3.93</td><td align="right">-1.82</td><td align="right">20.00</td><td align="right">10.00</td></tr>
<tr><td align="right">0.60</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.12</td><td align="right">4.31</td><td align="right">-2.19</td><td align="right">15.42</td><td align="right">9.86</td></tr>
<tr><td align="right">0.70</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.98</td><td align="right">4.67</td><td align="right">-2.69</td><td align="right">21.25</td><td align="right">8.89</td></tr>
<tr><td align="right">0.80</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.03</td><td align="right">4.87</td><td align="right">-2.84</td><td align="right">17.08</td><td align="right">8.33</td></tr>
<tr><td align="right">0.90</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">2.10</td><td align="right">5.15</td><td align="right">-3.05</td><td align="right">17.50</td><td align="right">8.47</td></tr>
<tr><td align="right">1.00</td><td align="right">8</td><td align="right">8</td><td align="right">0</td><td align="right">0</td><td align="right">0.007812</td><td align="right">1.89</td><td align="right">5.36</td><td align="right">-3.48</td><td align="right">20.00</td><td align="right">8.89</td></tr>
</tbody>
</table>

