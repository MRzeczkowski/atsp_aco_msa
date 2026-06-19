package workers

import (
	"atsp_aco_msa/modules/project"
	"os"
	"sync/atomic"
	"testing"
	"time"
)

func TestRunBoundedInstanceJobsRespectsWorkerLimit(t *testing.T) {
	atspsData := []project.AtspData{
		project.MakeAtspDataInResultsDirectory("a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		project.MakeAtspDataInResultsDirectory("b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		project.MakeAtspDataInResultsDirectory("c.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		project.MakeAtspDataInResultsDirectory("d.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		project.MakeAtspDataInResultsDirectory("e.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
	}

	var activeWorkers int32
	var maxActiveWorkers int32
	var processedJobs int32

	err := RunBoundedInstanceJobs(atspsData, 2, func(atspData project.AtspData) error {
		currentActiveWorkers := atomic.AddInt32(&activeWorkers, 1)
		for {
			currentMax := atomic.LoadInt32(&maxActiveWorkers)
			if currentActiveWorkers <= currentMax || atomic.CompareAndSwapInt32(&maxActiveWorkers, currentMax, currentActiveWorkers) {
				break
			}
		}

		time.Sleep(10 * time.Millisecond)
		atomic.AddInt32(&processedJobs, 1)
		atomic.AddInt32(&activeWorkers, -1)
		return nil
	})
	if err != nil {
		t.Fatalf("RunBoundedInstanceJobs returned unexpected error: %v", err)
	}

	if processedJobs != int32(len(atspsData)) {
		t.Fatalf("expected %d processed jobs, got %d", len(atspsData), processedJobs)
	}
	if maxActiveWorkers > 2 {
		t.Fatalf("expected at most two active workers, got %d", maxActiveWorkers)
	}
}

func TestRunBoundedIndexJobsRespectsWorkerLimitAndKeepsIndexedResults(t *testing.T) {
	result := make([]int, 6)
	var activeWorkers int32
	var maxActiveWorkers int32

	err := RunBoundedIndexJobs(len(result), 3, func(index int) error {
		currentActiveWorkers := atomic.AddInt32(&activeWorkers, 1)
		for {
			currentMax := atomic.LoadInt32(&maxActiveWorkers)
			if currentActiveWorkers <= currentMax || atomic.CompareAndSwapInt32(&maxActiveWorkers, currentMax, currentActiveWorkers) {
				break
			}
		}

		time.Sleep(10 * time.Millisecond)
		result[index] = index + 100
		atomic.AddInt32(&activeWorkers, -1)
		return nil
	})
	if err != nil {
		t.Fatalf("RunBoundedIndexJobs returned unexpected error: %v", err)
	}

	if maxActiveWorkers > 3 {
		t.Fatalf("expected at most three active workers, got %d", maxActiveWorkers)
	}
	for index, value := range result {
		if value != index+100 {
			t.Fatalf("expected indexed result %d to be %d, got %d", index, index+100, value)
		}
	}
}

func TestRunBoundedIndexJobsReturnsFirstError(t *testing.T) {
	err := RunBoundedIndexJobs(3, 1, func(index int) error {
		if index == 1 {
			return os.ErrInvalid
		}
		return nil
	})
	if err != os.ErrInvalid {
		t.Fatalf("expected %v, got %v", os.ErrInvalid, err)
	}
}

func TestRunIndexJobsWithSharedWorkersRespectsLimit(t *testing.T) {
	workerGate := make(chan struct{}, 2)
	var activeWorkers int32
	var maxActiveWorkers int32

	err := RunIndexJobsWithSharedWorkers(8, workerGate, func(index int) error {
		active := atomic.AddInt32(&activeWorkers, 1)
		for {
			currentMax := atomic.LoadInt32(&maxActiveWorkers)
			if active <= currentMax || atomic.CompareAndSwapInt32(&maxActiveWorkers, currentMax, active) {
				break
			}
		}

		time.Sleep(10 * time.Millisecond)
		atomic.AddInt32(&activeWorkers, -1)
		return nil
	})
	if err != nil {
		t.Fatalf("RunIndexJobsWithSharedWorkers returned unexpected error: %v", err)
	}
	if maxActiveWorkers > 2 {
		t.Fatalf("expected at most two active shared workers, got %d", maxActiveWorkers)
	}
}

func TestRunIndexJobsWithSharedWorkersReturnsFirstError(t *testing.T) {
	err := RunIndexJobsWithSharedWorkers(3, make(chan struct{}, 1), func(index int) error {
		if index == 1 {
			return os.ErrInvalid
		}
		return nil
	})
	if err != os.ErrInvalid {
		t.Fatalf("expected %v, got %v", os.ErrInvalid, err)
	}
}

func TestRunBoundedInstanceJobsReturnsFirstError(t *testing.T) {
	atspsData := []project.AtspData{
		project.MakeAtspDataInResultsDirectory("ok.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		project.MakeAtspDataInResultsDirectory("bad.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
	}

	err := RunBoundedInstanceJobs(atspsData, 1, func(atspData project.AtspData) error {
		if atspData.Name == "bad" {
			return os.ErrInvalid
		}
		return nil
	})
	if err != os.ErrInvalid {
		t.Fatalf("expected %v, got %v", os.ErrInvalid, err)
	}
}

func TestRunBoundedInstanceJobsHandlesEmptyInput(t *testing.T) {
	called := false
	err := RunBoundedInstanceJobs(nil, 2, func(atspData project.AtspData) error {
		called = true
		return nil
	})
	if err != nil {
		t.Fatalf("RunBoundedInstanceJobs returned unexpected error: %v", err)
	}
	if called {
		t.Fatal("job should not be called for empty input")
	}
}
