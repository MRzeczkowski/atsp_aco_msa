# Porównanie wyników GKS

Celem porównania było sprawdzenie, czy zaimplementowany wariant GKS daje wyniki zbliżone do tych z tabeli 3 w pracy:

Fred Glover, Gregory Gutin, Anders Yeo, Alexey Zverovich, "Construction heuristics for the asymmetric TSP", European Journal of Operational Research 129 (2001), 555-568.

W tabeli porównano kolumnę `GKS (%)` z pracy z wariantem `MSA patch bias = 0.00` z pliku `final/gks_deviation.md`. Ten wariant nie korzysta z informacji z MSA, więc odpowiada czystemu wariantowi GKS opartemu wyłącznie na koszcie łatania.

W pracy podano 26 instancji ATSP z TSPLIB. W lokalnym zestawie brakuje instancji `kro124p`, dlatego porównanie obejmuje 25 wspólnych instancji. Dla 12 z nich wyniki są zgodne z dokładnością do 0.01 p.p., dla 6 różnice są małe, a dla 7 wyraźne.

## Wyniki

| Instancja | GKS w pracy [%] | GKS w projekcie [%] | Różnica [p.p.] |
|---|---:|---:|---:|
| br17 | 0.00 | 0.00 | +0.00 |
| p43 | 0.32 | 0.27 | -0.05 |
| ry48p | 4.52 | 4.52 | +0.00 |
| ft53 | 12.31 | 12.31 | +0.00 |
| ft70 | 2.84 | 2.84 | +0.00 |
| ftv33 | 8.09 | 6.69 | -1.40 |
| ftv35 | 1.09 | 1.15 | +0.06 |
| ftv38 | 1.05 | 1.11 | +0.06 |
| ftv44 | 5.33 | 5.33 | +0.00 |
| ftv47 | 1.69 | 1.69 | +0.00 |
| ftv55 | 3.05 | 3.92 | +0.87 |
| ftv64 | 2.61 | 2.56 | -0.05 |
| ftv70 | 2.87 | 2.92 | +0.05 |
| ftv100 | 5.31 | 5.31 | +0.00 |
| ftv110 | 5.67 | 4.24 | -1.43 |
| ftv120 | 5.12 | 3.83 | -1.29 |
| ftv130 | 4.90 | 3.68 | -1.22 |
| ftv140 | 4.67 | 3.51 | -1.16 |
| ftv150 | 4.33 | 3.26 | -1.07 |
| ftv160 | 1.49 | 1.49 | +0.00 |
| ftv170 | 1.38 | 1.34 | -0.04 |
| kro124p | 8.69 | N/A | N/A |
| rbg323 | 0.00 | 0.00 | +0.00 |
| rbg358 | 0.00 | 0.00 | +0.00 |
| rbg403 | 0.00 | 0.00 | +0.00 |
| rbg443 | 0.00 | 0.00 | +0.00 |

## Średnie

| Zestaw | Średnie odchylenie [%] |
|---|---:|
| Praca, 26 instancji z tabeli 3 | 3.36 |
| Praca, 25 wspólnych instancji | 3.15 |
| Projekt, 25 wspólnych instancji | 2.88 |
| Projekt, wszystkie znane instancje | 6.57 |

Średnia dla wszystkich znanych instancji z projektu nie jest bezpośrednio porównywalna z wynikiem z pracy, ponieważ obejmuje dodatkowe instancje spoza tabeli 3.

## Uwagi

Wyniki nie są identyczne, ale duża część z nich się pokrywa. Różnice mogą wynikać z rozstrzygania remisów przy wyznaczaniu minimalnego pokrycia cyklowego oraz przy wyborze łatek. Praca nie podaje szczegółów pozwalających jednoznacznie odtworzyć te decyzje.
