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
1. Algorytm Edmonds'a na wersję Gabowa? 

2. Porozmawiać z Panem Solarzem o wykorzytaniu komputerów z uczelni? X

3. Analiza stytystyczna wpływu parametrów
4. Cykl Hamiltonowski w CMSA

5. Pan Krzysztof Sęp? X


6. Feromony niezależne od aktualnej heurystyki OK
7. Może wykres zbieżności? OK
8. 3-opt i sąsiedzi drugiego rzędu? OK

9. Latex vs Google docs

Ważne rzeczy w pracy:
1. TSP
2. MST
3. ATSP
4. MSA
5. ACO/MMAS
6. reduced 3-opt
