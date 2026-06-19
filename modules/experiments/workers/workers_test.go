package workers

import (
	"os"
	"sync/atomic"
	"testing"
	"time"
)

func TestRunJobsRespectsWorkerLimitAndKeepsIndexedResults(t *testing.T) {
	results := make([]int, 6)
	jobs := make([]Job, len(results))
	var activeWorkers int32
	var maxActiveWorkers int32

	for index := range results {
		index := index
		jobs[index] = Job{
			Label: "test",
			Run: func() error {
				currentActiveWorkers := atomic.AddInt32(&activeWorkers, 1)
				for {
					currentMax := atomic.LoadInt32(&maxActiveWorkers)
					if currentActiveWorkers <= currentMax || atomic.CompareAndSwapInt32(&maxActiveWorkers, currentMax, currentActiveWorkers) {
						break
					}
				}

				time.Sleep(10 * time.Millisecond)
				results[index] = index + 100
				atomic.AddInt32(&activeWorkers, -1)
				return nil
			},
		}
	}

	if err := RunJobs(jobs, 3); err != nil {
		t.Fatalf("RunJobs returned unexpected error: %v", err)
	}

	if maxActiveWorkers > 3 {
		t.Fatalf("expected at most three active workers, got %d", maxActiveWorkers)
	}
	for index, value := range results {
		if value != index+100 {
			t.Fatalf("expected indexed result %d to be %d, got %d", index, index+100, value)
		}
	}
}

func TestRunJobsReturnsFirstError(t *testing.T) {
	jobs := []Job{
		{Label: "ok", Run: func() error { return nil }},
		{Label: "bad", Run: func() error { return os.ErrInvalid }},
		{Label: "skipped", Run: func() error { return nil }},
	}

	if err := RunJobs(jobs, 1); err != os.ErrInvalid {
		t.Fatalf("expected %v, got %v", os.ErrInvalid, err)
	}
}

func TestRunJobsHandlesEmptyInput(t *testing.T) {
	if err := RunJobs(nil, 2); err != nil {
		t.Fatalf("RunJobs returned unexpected error: %v", err)
	}
}

func TestRunJobsRejectsInvalidWorkerCount(t *testing.T) {
	err := RunJobs([]Job{{Label: "test", Run: func() error { return nil }}}, 0)
	if err == nil {
		t.Fatal("expected invalid worker count error")
	}
}
