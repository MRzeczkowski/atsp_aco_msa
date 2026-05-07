# CMSA-Solution Analysis

Main high-CMSA threshold: 0.8

## Summary

- Instances analyzed: 28
- Instances with found optimal tours: 23
- Instances without found optimal tours: 5
- Average optimal-edge CMSA coverage: 60.54%
- Average high-CMSA precision: 58.40%
- Average high-CMSA recall: 40.65%
- Average high-CMSA lift: 43.505372

## Notes

- `solutions.csv` is treated as the set of found optimal tours; `Opt edges` is the directed edge union across those tours.
- `Opt edge density` is the random edge-hit baseline used by lift.
- Normalized CMSA values use `CMSA / (dimension - 1)`, matching the scale used by the heuristic.
- Average and median normalized CMSA on found-optimal edges include missing edges as zero.
- Rows with no found optimal tours keep CMSA density metrics, but overlap metrics are set to zero because there is no found-tour edge set.
- Edges marked as not seen in found optimal tours are not proven non-optimal; they are absent from the tours currently recorded in `solutions.csv`.
- `Lift` is high-CMSA precision divided by `Opt edge density`, so it can be large when the found-optimal edge union is sparse.

## Instances

| Instance | Unique tours | Opt edges | Opt edge density % | CMSA density % | Coverage % | High precision % | High recall % | Lift | Missing | High not seen |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| atex5 | 158 | 353 | 6.91 | 2.72 | 34.84 | 88.06 | 16.71 | 12.752442 | 230 | 8 |
| code198 | 1482 | 27073 | 69.41 | 1.01 | 0.30 | 36.55 | 0.27 | 0.526576 | 26991 | 125 |
| crane100_0 | 1 | 100 | 1.01 | 1.77 | 66.00 | 44.21 | 42.00 | 43.768421 | 34 | 53 |
| crane100_1 | 1 | 100 | 1.01 | 1.77 | 68.00 | 53.33 | 48.00 | 52.800000 | 32 | 42 |
| crane100_2 | 1 | 100 | 1.01 | 1.68 | 68.00 | 54.65 | 47.00 | 54.104651 | 32 | 39 |
| crane66_0 | 1 | 66 | 1.54 | 2.66 | 69.70 | 51.72 | 45.45 | 33.620690 | 20 | 28 |
| crane66_1 | 1 | 66 | 1.54 | 2.75 | 62.12 | 40.98 | 37.88 | 26.639344 | 25 | 36 |
| crane66_2 | 1 | 66 | 1.54 | 2.66 | 74.24 | 61.67 | 56.06 | 40.083333 | 17 | 23 |
| dc112 | 0 | 0 | 0.00 | 1.78 | 0.00 | 0.00 | 0.00 | 0.000000 | 0 | 0 |
| dc126 | 0 | 0 | 0.00 | 1.61 | 0.00 | 0.00 | 0.00 | 0.000000 | 0 | 0 |
| dc134 | 0 | 0 | 0.00 | 1.78 | 0.00 | 0.00 | 0.00 | 0.000000 | 0 | 0 |
| dc176 | 0 | 0 | 0.00 | 1.28 | 0.00 | 0.00 | 0.00 | 0.000000 | 0 | 0 |
| dc188 | 0 | 0 | 0.00 | 1.06 | 0.00 | 0.00 | 0.00 | 0.000000 | 0 | 0 |
| ft53 | 1 | 53 | 1.92 | 3.23 | 69.81 | 45.10 | 43.40 | 23.450980 | 16 | 28 |
| ft70 | 1 | 70 | 1.45 | 2.42 | 58.57 | 48.48 | 45.71 | 33.454545 | 29 | 34 |
| ftv100 | 81 | 135 | 1.34 | 1.54 | 65.19 | 60.82 | 43.70 | 45.505918 | 47 | 38 |
| ftv110 | 63 | 142 | 1.16 | 1.40 | 66.20 | 61.90 | 45.77 | 53.229376 | 48 | 40 |
| ftv120 | 64 | 158 | 1.09 | 1.27 | 63.92 | 61.74 | 44.94 | 56.737479 | 57 | 44 |
| ftv130 | 133 | 173 | 1.02 | 1.15 | 65.90 | 65.60 | 47.40 | 64.576185 | 59 | 43 |
| ftv140 | 109 | 183 | 0.93 | 1.09 | 67.21 | 65.19 | 48.09 | 70.314511 | 60 | 47 |
| ftv150 | 148 | 201 | 0.89 | 1.01 | 66.17 | 69.44 | 49.75 | 78.254561 | 68 | 44 |
| ftv160 | 149 | 218 | 0.85 | 0.95 | 67.43 | 68.92 | 46.79 | 81.438135 | 71 | 46 |
| ftv170 | 531 | 232 | 0.80 | 0.88 | 65.95 | 65.61 | 44.40 | 82.204316 | 79 | 54 |
| ftv55 | 111 | 80 | 2.60 | 2.82 | 65.00 | 61.70 | 36.25 | 23.755319 | 28 | 18 |
| ftv64 | 119 | 92 | 2.21 | 2.38 | 64.13 | 62.07 | 39.13 | 28.065967 | 33 | 22 |
| ftv70 | 107 | 93 | 1.87 | 2.15 | 67.74 | 59.38 | 40.86 | 31.730511 | 30 | 26 |
| ftv90 | 40 | 113 | 1.38 | 1.72 | 68.14 | 59.77 | 46.02 | 43.320110 | 36 | 35 |
| td100_1 | 100 | 280 | 2.77 | 1.74 | 27.86 | 56.25 | 19.29 | 20.290179 | 202 | 42 |
