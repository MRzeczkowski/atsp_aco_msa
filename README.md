# ATSP ACO MSA

Program uruchamia eksperymenty dla asymetrycznego problemu komiwojażera. Porównuje wariant bazowy MMAS z heurystykami opartymi na MSA i pokryciu cyklowym.

## Wymagania

Potrzebne są:

- Go w wersji 1.22.1 lub nowszej;
- terminal;
- pliki instancji ATSP w katalogu `tsplib_files`.

Wymagane instancje są już zapisane w repozytorium. Pozostałe biblioteki pobiera narzędzie Go.

## Instalacja Go

1. Pobierz instalator ze strony <https://go.dev/dl/>.
2. Zainstaluj Go.
3. Otwórz terminal.
4. Sprawdź instalację:

```bash
go version
```

Polecenie powinno wyświetlić wersję Go. Jeżeli system nie rozpoznaje polecenia `go`, uruchom ponownie terminal albo dodaj Go do zmiennej `PATH`.

## Przygotowanie projektu

W terminalu przejdź do głównego katalogu repozytorium. Musi to być katalog zawierający pliki `go.mod` i `main.go`.

Pobierz zależności:

```bash
go mod download
```

Uruchom testy:

```bash
go test ./...
```

## Kompilacja

macOS lub Linux:

```bash
go build -o atsp_aco_msa .
```

Windows:

```powershell
go build -o atsp_aco_msa.exe .
```

Po kompilacji plik wykonywalny znajduje się w głównym katalogu repozytorium.

## Uruchamianie

Program należy uruchamiać z głównego katalogu repozytorium. Korzysta on ze ścieżek względnych do katalogów `tsplib_files` i `artifacts`.

Bez kompilacji:

```bash
go run . [flagi]
```

Po kompilacji na macOS lub Linux:

```bash
./atsp_aco_msa [flagi]
```

Po kompilacji na Windows:

```powershell
.\atsp_aco_msa.exe [flagi]
```

Lista wszystkich flag:

```bash
go run . -h
```

Samo `go run .` uruchamia strojenie wszystkich heurystyk na całym zbiorze `tuning`. Do pierwszego sprawdzenia programu użyj krótszego polecenia z następnej sekcji.

## Szybki test

To polecenie uruchamia strojenie jednej heurystyki dla małej instancji `br17` przy użyciu jednego zadania wykonywanego w danym momencie:

```bash
go run . -mode experiment -instances smoke -heuristic strict-msa -workers 1
```

Program automatycznie utworzy brakujące dane MSA i zapisze wynik w katalogu `artifacts`.

## Tryby pracy

Tryb wybiera flaga `-mode`.

### `experiment`

Uruchamia strojenie parametrów heurystyk bez 3-opt. Dla każdego zestawu parametrów wykonuje 30 uruchomień MMAS.

```bash
go run . -mode experiment -instances tuning
```

Wyniki są zapisywane w:

```text
artifacts/experiments/tuning
```

Jeżeli flaga `-heuristic` zostanie pominięta, program sprawdzi wszystkie cztery heurystyki dostępne w tym trybie.

### `evaluation`

Uruchamia końcową ewaluację bez 3-opt. Dla każdej konfiguracji wykonuje 50 uruchomień MMAS. Domyślnie wykorzystuje zbiór instancji `evaluation`.

```bash
go run . -mode evaluation
```

Wyniki są zapisywane w:

```text
artifacts/experiments/evaluation/no_3opt
```

Po wykonaniu wszystkich podstawowych wariantów program automatycznie uruchamia analizę wyników.

### `evaluation+3opt`

Działa tak samo jak `evaluation`, ale po zbudowaniu cyklu stosuje zredukowany 3-opt.

```bash
go run . -mode evaluation+3opt
```

Wyniki są zapisywane w:

```text
artifacts/experiments/evaluation/with_3opt
```

### `analyze`

Nie uruchamia eksperymentów MMAS. Odczytuje istniejące wyniki i tworzy raporty, podsumowania oraz wykresy.

Pełna analiza wyników ewaluacji:

```bash
go run . -mode analyze -instances evaluation -analysis all
```

Samo podsumowanie strojenia:

```bash
go run . -mode analyze -analysis tuning
```

Sama analiza odchylenia GKS:

```bash
go run . -mode analyze -analysis gks-deviation
```

### `all`

Najpierw uruchamia tryb `experiment`, a następnie `analyze` dla wybranego zbioru instancji. Nie uruchamia trybów `evaluation` ani `evaluation+3opt`.

```bash
go run . -mode all -instances tuning
```

### `rebuild-cache`

Ponownie wyznacza MSA i minimalne pokrycia cyklowe, a następnie odtwarza związane z nimi wykresy. Bez flagi `-instances` przetwarza wszystkie instancje ze znanym optimum.

```bash
go run . -mode rebuild-cache
```

Ten tryb nadpisuje dane w katalogach:

```text
artifacts/cache/msa
artifacts/cache/cycle_cover
```

## Zbiory instancji

Zbiór wybiera flaga `-instances`:

- `smoke` - tylko mała instancja `br17`, przeznaczona do szybkiego sprawdzenia programu;
- `tuning` - instancje używane do strojenia parametrów;
- `evaluation` - instancje używane do końcowej ewaluacji;
- `all-known` - wszystkie dostępne instancje ze znanym optimum.

Przykład:

```bash
go run . -mode rebuild-cache -instances smoke
```

## Wybór heurystyki podczas strojenia

Flaga `-heuristic` działa w trybach `experiment` i `all`.

Dostępne wartości:

- `all` - wszystkie poniższe heurystyki;
- `strict-msa` - Heurystyka MSA;
- `rooted-msa` - Zakorzenione MSA;
- `cycle-cover` - heurystyka pokrycia cyklowego;
- `cycle-cover-msa-patching` - GKS+MSA.

Przykład:

```bash
go run . -mode experiment -instances tuning -heuristic rooted-msa
```

Pominięcie flagi oznacza to samo co `-heuristic all`.

## Wybór wariantu podczas ewaluacji

Flaga `-evaluation-heuristic` działa w trybach `evaluation` i `evaluation+3opt`.

Dostępne wartości:

- `all` - wszystkie podstawowe warianty;
- `controls` - wszystkie warianty kontrolne;
- `baseline`;
- `strict-msa`;
- `rooted-msa`;
- `random-sparse`;
- `distance-ranked-sparse`;
- `shuffled-msa`;
- `cycle-cover`;
- `cycle-cover-patching`;
- `cycle-cover-msa-patching`.

Przykład uruchomienia samych kontroli:

```bash
go run . -mode evaluation -evaluation-heuristic controls
```

## Liczba równoległych zadań

Flaga `-workers` określa maksymalną liczbę konfiguracji wykonywanych jednocześnie.

- `-workers 1` - wykonywanie sekwencyjne;
- `-workers 4` - maksymalnie cztery zadania jednocześnie;
- `-workers 0` albo brak flagi - połowa dostępnych procesorów logicznych.

Przykład:

```bash
go run . -mode evaluation -workers 4
```

Większa wartość może skrócić czas obliczeń, ale zwiększa wykorzystanie procesora i pamięci.

## Pliki wynikowe

Program zapisuje wszystkie wygenerowane dane w katalogu `artifacts`:

- `artifacts/cache/msa` - drzewa MSA, macierze heurystyki i wykresy;
- `artifacts/cache/cycle_cover` - minimalne pokrycia cyklowe i wykresy;
- `artifacts/experiments/tuning` - wyniki strojenia;
- `artifacts/experiments/evaluation/no_3opt` - ewaluacja bez 3-opt;
- `artifacts/experiments/evaluation/with_3opt` - ewaluacja z 3-opt;
- `artifacts/experiments/evaluation/controls` - wyniki wariantów kontrolnych;
- `artifacts/solutions` - znalezione cykle optymalne i ich analiza.

Ponowne uruchomienie eksperymentu może nadpisać istniejące wyniki dla tej samej konfiguracji.
