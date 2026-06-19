package workers

import (
	"fmt"
	"sync"
)

type Job struct {
	Label string
	Run   func() error
}

func RunJobs(jobs []Job, workers int) error {
	if len(jobs) == 0 {
		return nil
	}
	if workers < 1 {
		return fmt.Errorf("workers must be at least one")
	}

	effectiveWorkers := min(workers, len(jobs))
	nextIndex := 0
	var firstErr error
	var mutex sync.Mutex
	var waitGroup sync.WaitGroup

	runWorker := func() {
		defer waitGroup.Done()

		for {
			mutex.Lock()
			if firstErr != nil || nextIndex >= len(jobs) {
				mutex.Unlock()
				return
			}
			index := nextIndex
			nextIndex++
			mutex.Unlock()

			if err := jobs[index].Run(); err != nil {
				mutex.Lock()
				if firstErr == nil {
					firstErr = err
				}
				mutex.Unlock()
				return
			}
		}
	}

	waitGroup.Add(effectiveWorkers)
	for i := 0; i < effectiveWorkers; i++ {
		go runWorker()
	}
	waitGroup.Wait()

	return firstErr
}
