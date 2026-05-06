# Cycle-Cover/CMSA Analysis

Main high-CMSA threshold: 1.0

## Summary

- Instances analyzed: 28
- Instances with found optimal tours: 23
- Average cycle-cover precision: 73.77%
- Average cycle-cover recall: 53.17%
- Average high-CMSA precision: 67.67%
- Average high-CMSA recall: 29.98%
- Average cycle-cover-gated high-CMSA precision: 84.48%
- Average cycle-cover-gated high-CMSA recall: 24.91%
- Average precision gain from cycle-cover gate: 0.168140
- Average recall loss from cycle-cover gate: 0.050730

## Notes

- `Cycle cover` is the minimum assignment/cycle-cover solution computed from the original ATSP matrix with self-loops forbidden.
- `High CMSA` uses normalized `CMSA / (dimension - 1)` and the threshold shown above.
- `Cycle-cover high-CMSA` is the strict intersection of high-CMSA edges and minimum-cycle-cover edges.
- Precision/recall use the union of tours recorded in `solutions.csv`; absent edges are not proven non-optimal.
- Precision gain and recall loss compare `cycle-cover high-CMSA` against high-CMSA alone.

## Instances

| Instance | Tours | Opt edges | CC gap % | CC cycles | CC precision % | CC recall % | High precision % | High recall % | Gated precision % | Gated recall % | Precision gain | Recall loss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| atex5 | 156 | 353 | 97.57 | 35 | 100.00 | 20.40 | 100.00 | 9.07 | 100.00 | 9.07 | 0.000000 | 0.000000 |
| code198 | 915 | 26994 | 0.09 | 4 | 95.45 | 0.70 | 37.97 | 0.26 | 50.00 | 0.00 | 0.120321 | 0.002593 |
| crane100_0 | 1 | 100 | 8.55 | 22 | 53.00 | 53.00 | 53.57 | 30.00 | 64.10 | 25.00 | 0.105311 | 0.050000 |
| crane100_1 | 1 | 100 | 5.57 | 17 | 60.00 | 60.00 | 57.89 | 33.00 | 84.38 | 27.00 | 0.264803 | 0.060000 |
| crane100_2 | 1 | 100 | 10.30 | 29 | 58.00 | 58.00 | 63.33 | 38.00 | 82.86 | 29.00 | 0.195238 | 0.090000 |
| crane66_0 | 1 | 66 | 12.13 | 15 | 40.91 | 40.91 | 61.11 | 33.33 | 66.67 | 21.21 | 0.055556 | 0.121212 |
| crane66_1 | 1 | 66 | 16.84 | 22 | 54.55 | 54.55 | 54.55 | 27.27 | 80.00 | 24.24 | 0.254545 | 0.030303 |
| crane66_2 | 1 | 66 | 9.55 | 22 | 48.48 | 48.48 | 73.68 | 42.42 | 82.61 | 28.79 | 0.089245 | 0.136364 |
| dc112 | 0 | 0 | 0.89 | 25 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.000000 | 0.000000 |
| dc126 | 0 | 0 | 3.81 | 51 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.000000 | 0.000000 |
| dc134 | 0 | 0 | 0.37 | 21 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.000000 | 0.000000 |
| dc176 | 0 | 0 | 0.84 | 41 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.000000 | 0.000000 |
| dc188 | 0 | 0 | 1.24 | 30 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.000000 | 0.000000 |
| ft53 | 1 | 53 | 14.11 | 8 | 64.15 | 64.15 | 53.57 | 28.30 | 69.23 | 16.98 | 0.156593 | 0.113208 |
| ft70 | 1 | 70 | 1.80 | 10 | 65.71 | 65.71 | 51.92 | 38.57 | 96.30 | 37.14 | 0.443732 | 0.014286 |
| ftv100 | 81 | 135 | 6.60 | 12 | 86.14 | 64.44 | 74.14 | 31.85 | 90.70 | 28.89 | 0.165597 | 0.029630 |
| ftv110 | 62 | 142 | 5.87 | 13 | 79.28 | 61.97 | 68.75 | 30.99 | 84.44 | 26.76 | 0.156944 | 0.042254 |
| ftv120 | 64 | 158 | 6.51 | 16 | 76.86 | 58.86 | 68.49 | 31.65 | 86.54 | 28.48 | 0.180453 | 0.031646 |
| ftv130 | 133 | 173 | 4.46 | 16 | 80.92 | 61.27 | 75.00 | 34.68 | 91.53 | 31.21 | 0.165254 | 0.034682 |
| ftv140 | 109 | 183 | 4.92 | 17 | 80.85 | 62.30 | 75.29 | 34.97 | 90.48 | 31.15 | 0.151821 | 0.038251 |
| ftv150 | 148 | 201 | 3.91 | 16 | 84.11 | 63.18 | 76.09 | 34.83 | 92.42 | 30.35 | 0.163373 | 0.044776 |
| ftv160 | 149 | 218 | 4.58 | 16 | 88.82 | 65.60 | 76.60 | 33.03 | 92.86 | 29.82 | 0.162614 | 0.032110 |
| ftv170 | 531 | 232 | 4.50 | 20 | 82.46 | 60.78 | 76.24 | 33.19 | 90.54 | 28.88 | 0.143029 | 0.043103 |
| ftv55 | 111 | 80 | 10.76 | 9 | 73.21 | 51.25 | 80.00 | 30.00 | 86.96 | 25.00 | 0.069565 | 0.050000 |
| ftv64 | 119 | 92 | 6.42 | 13 | 73.85 | 52.17 | 76.32 | 31.52 | 92.59 | 27.17 | 0.162768 | 0.043478 |
| ftv70 | 107 | 93 | 9.44 | 9 | 70.42 | 53.76 | 71.43 | 32.26 | 78.12 | 26.88 | 0.066964 | 0.053763 |
| ftv90 | 40 | 113 | 6.33 | 11 | 82.42 | 66.37 | 73.08 | 33.63 | 89.74 | 30.97 | 0.166667 | 0.026549 |
| td100_1 | 100 | 280 | 0.00 | 8 | 97.03 | 35.00 | 57.32 | 16.79 | 100.00 | 8.93 | 0.426829 | 0.078571 |
