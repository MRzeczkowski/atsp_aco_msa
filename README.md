# Hybrid algorithm for ATSP

## ATSP

- Paper detailing ATSP: https://www.diva-portal.org/smash/get/diva2:358638/fulltext01.pdf

- Analysis of Heuristics for the ATSP: https://www.researchgate.net/publication/2404400_Experimental_Analysis_of_Heuristics_for_the_ATSP

- Asadpours paper on approximation: https://homes.cs.washington.edu/~shayan/atsp.pdf

- NetworkX on ATSP: https://blog.scientific-python.org/posts/networkx/atsp/

- TSPLib test data: http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/atsp/

## MSA

- https://en.wikipedia.org/wiki/Arborescence_(graph_theory)

- https://en.wikipedia.org/wiki/Edmonds%27_algorithm

- https://en.wikipedia.org/wiki/Multitree

- MSA presentation: https://homes.di.unimi.it/righini/Didattica/OttimizzazioneCombinatoria/MaterialeOC/5%20-%20Minimum%20spanning%20r-arborescence.pdf

- Example in Python with some explanations: https://wendy-xiao.github.io/posts/2020-07-10-chuliuemdond_algorithm/

- Gabow et. al. improved version of Edmonds' algorithm: https://sci-hub.se/10.1007/BF02579168

- Implementation with ideas from Trajan, maybe also Gabow: https://github.com/atofigh/edmonds-alg

- Disjoint-set implementation in Go: https://github.com/spakin/disjoint

## ACO

- Hybrid ACS for TSP: https://sci-hub.se/10.1109/ICICTA.2010.731

- MMAS paper: https://sci-hub.se/https://doi.org/10.1016/S0167-739X(00)00043-1

- Another paper on MMAS, with info on 3-opt: https://lia.disi.unibo.it/Courses/SistInt/articoli/max-min-ant.pdf

- ACO for TSP: https://faculty.washington.edu/paymana/swarm/stutzle99-eaecs.pdf

- ACO for TSP parameter analysis: http://article.nadiapub.com/IJUNESST/vol7_no4/16.pdf

- More on parameter tuning: https://sci-hub.se/https://doi.org/10.1109/CSNT.2018.8820263

- Something about pheromones: https://sci-hub.se/10.1109/PIC.2014.6972311

- Model induced MMAS: https://sci-hub.se/https://doi.org/10.1016/j.asoc.2012.04.008

- ACO variants overview: https://sci-hub.se/http://dx.doi.org/10.17485/ijst/2015/v8i31/87296

- Ant-Q: https://sci-hub.se/http://dx.doi.org/10.1016/b978-1-55860-377-6.50039-6

- Ant-Q presentation: https://csc.csudh.edu/btang/seminar/PDDC/Ant-Q.pdf

## 3-opt

- Algorithms for TSP, some info on 3-opt: https://sci-hub.se/https://doi.org/10.1287/ijoc.4.4.387

- Local optimizations for TSP, 3.3 is useful: https://www.cs.ubc.ca/~hutter/previous-earg/EmpAlgReadingGroup/TSP-JohMcg97.pdf

- Markov chains for TSP: https://content.wolfram.com/sites/13/2018/02/05-3-3.pdf

- Implementation in C: https://github.com/ozanyerli/tsp3opt

- Implementation in Python: https://github.com/BraveDistribution/pytsp/blob/master/pytsp/k_opt_tsp.py

- Article related to the Python implementation: https://github.com/BraveDistribution/pytsp/blob/master/pytsp/k_opt_tsp.py

- Best (still bad) explanation of reduced 3-opt, 3.3 is useful: https://www.cs.ubc.ca/labs/algorithms/Courses/CPSC532D-03/Resources/StuHoo99.pdf

- Pointer to the reduced 3-opt algorithm was found here: https://www.es.ele.tue.nl/cps/publications/pinxten2016mogtsp.pdf

- Blog about TSP and in big part about k-opt algorithms: https://github.com/BraveDistribution/pytsp/blob/master/pytsp/k_opt_tsp.py

## Others

- "Ghouila-Houri from 1960 asserts that every directed graph on n vertices with minimum out-degree and in-degree at least n/2 contains a directed Hamilton cycle."

- On Hamiltonian cycles in directed graphs: https://www.sciencedirect.com/science/article/pii/S0195669811001788


Co dalej?

6. Feromony niezależne od aktualnej heurystyki
 - niewiele dało, usunąłem, opisać po krótce w pracy

7. Może wykres zbieżności? 
 - OK: zrobione, ale bardzo dużo wykresów i raczej lepiej je wygenerować dla określonych ekperymentów, dodać przykład z jednego eksperymentu do uśrednionego wykresu

11. Zmniejszyć ilość mrówek?
 - losowe przypisanie mrówek do węzłów początkowych
 - ustalić ilość mrówek, raczej jako procent węzłów, sprawdzić zmienną ilość mrówek i zapisać numer najlepsze mrówki

12. Ile iteracji robi 3-opt?
 - dodać to wyników, zdaje się że CMSA skłania do globalnego optimum, ale po drodze do niego ścieżki są suboptymalne

13. FDC
 - policzyć!

Rozdziały:

1. TSP
    1. Opis problemu
    2. Związek z MST
        1. Algorytmy

2. ATSP
    1. Opis problemu - różnice względem TSP
    2. MSA
        1. Algorytm Edmondsa
    3. Niezbadany związek z MSA

3. ACO
    1. Zasada działania
    2. MMAS
    3. W połączeniu z innymi algorytmami (MSA + 3-opt)

4. Algorytm hybrydowy - mój
    1. CMSA
    2. Wykorzystanie CMSA
    3. "reduced 3-opt"

5. Metodyka
    1. Dane testowe
        1. FDC
    2. Eksperymenty

6. Wyniki

7. Podsumowanie


Normal:
Creating results/br17_cmsa.csv took: 10 ms
Creating results/ft53_cmsa.csv took: 384 ms
Creating results/ft70_cmsa.csv took: 1011 ms
Creating results/ftv170_cmsa.csv took: 51177 ms
Creating results/ftv33_cmsa.csv took: 100 ms
Creating results/ftv35_cmsa.csv took: 116 ms
Creating results/ftv38_cmsa.csv took: 145 ms
Creating results/ftv44_cmsa.csv took: 211 ms
Creating results/ftv47_cmsa.csv took: 313 ms
Creating results/ftv55_cmsa.csv took: 608 ms
Creating results/ftv64_cmsa.csv took: 1078 ms
Creating results/ftv70_cmsa.csv took: 1507 ms
Creating results/p43_cmsa.csv took: 623 ms
Creating results/rbg323_cmsa.csv took: 1764136 ms
Creating results/rbg358_cmsa.csv took: 2406854 ms
Creating results/rbg403_cmsa.csv took: 3645321 ms
Creating results/rbg443_cmsa.csv took: 5475710 ms
Creating results/ry48p_cmsa.csv took: 483 ms

Preallocation:
Creating results/br17_cmsa.csv took: 12 ms
Creating results/ft53_cmsa.csv took: 310 ms
Creating results/ft70_cmsa.csv took: 797 ms
Creating results/ftv170_cmsa.csv took: 40985 ms
Creating results/ftv33_cmsa.csv took: 79 ms
Creating results/ftv35_cmsa.csv took: 90 ms
Creating results/ftv38_cmsa.csv took: 112 ms
Creating results/ftv44_cmsa.csv took: 156 ms
Creating results/ftv47_cmsa.csv took: 227 ms
Creating results/ftv55_cmsa.csv took: 469 ms
Creating results/ftv64_cmsa.csv took: 815 ms
Creating results/ftv70_cmsa.csv took: 1115 ms
Creating results/p43_cmsa.csv took: 445 ms
Creating results/rbg323_cmsa.csv took: 1498863 ms
Creating results/rbg358_cmsa.csv took: 2061367 ms
Creating results/rbg403_cmsa.csv took: 3091203 ms
Creating results/rbg443_cmsa.csv took: 4772496 ms
Creating results/ry48p_cmsa.csv took: 369 ms

Grouping:
Creating results/br17_cmsa.csv took: 8 ms
Creating results/ft53_cmsa.csv took: 232 ms
Creating results/ft70_cmsa.csv took: 573 ms
Creating results/ftv170_cmsa.csv took: 26199 ms
Creating results/ftv33_cmsa.csv took: 62 ms
Creating results/ftv35_cmsa.csv took: 71 ms
Creating results/ftv38_cmsa.csv took: 95 ms
Creating results/ftv44_cmsa.csv took: 128 ms
Creating results/ftv47_cmsa.csv took: 184 ms
Creating results/ftv55_cmsa.csv took: 359 ms
Creating results/ftv64_cmsa.csv took: 645 ms
Creating results/ftv70_cmsa.csv took: 837 ms
Creating results/p43_cmsa.csv took: 353 ms
Creating results/rbg323_cmsa.csv took: 824796 ms
Creating results/rbg358_cmsa.csv took: 1085652 ms
Creating results/rbg403_cmsa.csv took: 1624324 ms
Creating results/rbg443_cmsa.csv took: 2496786 ms
Creating results/ry48p_cmsa.csv took: 305 ms