# MSA Support / Solution Analysis

- Instances analyzed: 28
- Instances with found optimal tours: 23
- Instances without found optimal tours: 5

## How To Read This

- `Precision` is the share of selected MSA support edges that were seen in found optimal tours.
- `Recall` is the share of the found-optimal edge union selected by the MSA support variant.
- `Lift` is precision divided by random edge-hit probability for that instance; higher means the selected edges are less random.
- `Tour coverage` is per-tour, not union-based: it measures how much of one complete found optimal tour is covered by the selected edges.

## Raw MSA Support Coverage

This checks whether MSA support contains found-optimal edges at all before selecting or filtering them.

| Metric | Average |
|---|---:|
| Found-optimal edges present in MSA support | 60.54% |

| Instance | Coverage % |
|---|---:|
| atex5 | 34.84 |
| code198 | 0.30 |
| crane100_0 | 66.00 |
| crane100_1 | 68.00 |
| crane100_2 | 68.00 |
| crane66_0 | 69.70 |
| crane66_1 | 62.12 |
| crane66_2 | 74.24 |
| ft53 | 69.81 |
| ft70 | 58.57 |
| ftv100 | 65.19 |
| ftv110 | 66.20 |
| ftv120 | 63.92 |
| ftv130 | 65.90 |
| ftv140 | 67.21 |
| ftv150 | 66.17 |
| ftv160 | 67.43 |
| ftv170 | 65.95 |
| ftv55 | 65.00 |
| ftv64 | 64.13 |
| ftv70 | 67.74 |
| ftv90 | 68.14 |
| td100_1 | 27.86 |

## Threshold MSA support

This variant selects every edge with normalized MSA support at least 0.8.

| Metric | Average |
|---|---:|
| Selected edges | 96.0 |
| Precision | 58.40% |
| Recall | 40.65% |
| Lift | 43.505272 |

| Instance | Edges | Precision % | Recall % | Lift |
|---|---:|---:|---:|---:|
| atex5 | 67 | 88.06 | 16.71 | 12.752442 |
| code198 | 197 | 36.55 | 0.26 | 0.524291 |
| crane100_0 | 95 | 44.21 | 42.00 | 43.768421 |
| crane100_1 | 90 | 53.33 | 48.00 | 52.800000 |
| crane100_2 | 86 | 54.65 | 47.00 | 54.104651 |
| crane66_0 | 58 | 51.72 | 45.45 | 33.620690 |
| crane66_1 | 61 | 40.98 | 37.88 | 26.639344 |
| crane66_2 | 60 | 61.67 | 56.06 | 40.083333 |
| ft53 | 51 | 45.10 | 43.40 | 23.450980 |
| ft70 | 66 | 48.48 | 45.71 | 33.454545 |
| ftv100 | 97 | 60.82 | 43.70 | 45.505918 |
| ftv110 | 105 | 61.90 | 45.77 | 53.229376 |
| ftv120 | 115 | 61.74 | 44.94 | 56.737479 |
| ftv130 | 125 | 65.60 | 47.40 | 64.576185 |
| ftv140 | 135 | 65.19 | 48.09 | 70.314511 |
| ftv150 | 144 | 69.44 | 49.75 | 78.254561 |
| ftv160 | 148 | 68.92 | 46.79 | 81.438135 |
| ftv170 | 157 | 65.61 | 44.40 | 82.204316 |
| ftv55 | 47 | 61.70 | 36.25 | 23.755319 |
| ftv64 | 58 | 62.07 | 39.13 | 28.065967 |
| ftv70 | 64 | 59.38 | 40.86 | 31.730511 |
| ftv90 | 87 | 59.77 | 46.02 | 43.320110 |
| td100_1 | 96 | 56.25 | 19.29 | 20.290179 |

## Top N-1 MSA support

This variant selects a fixed number of edges: the top `dimension - 1` positive MSA support edges, with deterministic tie-breaking by edge id.

| Metric | Average |
|---|---:|
| Selected edges | 101.7 |
| Precision | 56.51% |
| Recall | 41.95% |
| Lift | 42.195284 |
| Average tour coverage | 46.82% |
| Best-tour coverage | 48.62% |

| Instance | Edges | Precision % | Recall % | Avg tour coverage % | Best-tour coverage % |
|---|---:|---:|---:|---:|---:|
| atex5 | 71 | 85.92 | 17.28 | 37.23 | 44.44 |
| code198 | 197 | 36.55 | 0.26 | 1.01 | 1.01 |
| crane100_0 | 99 | 43.43 | 43.00 | 43.00 | 43.00 |
| crane100_1 | 99 | 50.51 | 50.00 | 50.00 | 50.00 |
| crane100_2 | 99 | 53.54 | 53.00 | 53.00 | 53.00 |
| crane66_0 | 65 | 46.15 | 45.45 | 45.45 | 45.45 |
| crane66_1 | 65 | 38.46 | 37.88 | 37.88 | 37.88 |
| crane66_2 | 65 | 56.92 | 56.06 | 56.06 | 56.06 |
| ft53 | 52 | 44.23 | 43.40 | 43.40 | 43.40 |
| ft70 | 69 | 46.38 | 45.71 | 45.71 | 45.71 |
| ftv100 | 100 | 60.00 | 44.44 | 50.40 | 53.47 |
| ftv110 | 110 | 61.82 | 47.89 | 52.84 | 54.95 |
| ftv120 | 120 | 61.67 | 46.84 | 53.76 | 55.37 |
| ftv130 | 130 | 64.62 | 48.55 | 56.08 | 58.02 |
| ftv140 | 140 | 64.29 | 49.18 | 56.77 | 58.16 |
| ftv150 | 150 | 68.00 | 50.75 | 59.94 | 62.25 |
| ftv160 | 160 | 66.88 | 49.08 | 58.37 | 60.25 |
| ftv170 | 170 | 63.53 | 46.55 | 54.51 | 57.31 |
| ftv55 | 55 | 56.36 | 38.75 | 45.46 | 50.00 |
| ftv64 | 64 | 59.38 | 41.30 | 49.54 | 53.85 |
| ftv70 | 70 | 57.14 | 43.01 | 48.46 | 52.11 |
| ftv90 | 90 | 58.89 | 46.90 | 50.93 | 53.85 |
| td100_1 | 100 | 55.00 | 19.64 | 27.07 | 28.71 |

Instances without found optimal tours are omitted from the overlap tables because precision, recall, and tour coverage cannot be interpreted without a found-tour edge set.
