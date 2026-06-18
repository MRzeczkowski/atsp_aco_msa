package app

import (
	"fmt"
	"sync"
)

func runBoundedInstanceJobs(atspsData []AtspData, workers int, job func(AtspData) error) error {
	if len(atspsData) == 0 {
		return nil
	}

	return runBoundedIndexJobs(len(atspsData), workers, func(index int) error {
		return job(atspsData[index])
	})
}

func runBoundedIndexJobs(jobCount, workers int, job func(int) error) error {
	if jobCount == 0 {
		return nil
	}
	if workers < 1 {
		return fmt.Errorf("workers must be at least one")
	}

	effectiveWorkers := min(workers, jobCount)
	nextIndex := 0
	var firstErr error
	var mutex sync.Mutex
	var waitGroup sync.WaitGroup

	runWorker := func() {
		defer waitGroup.Done()

		for {
			mutex.Lock()
			if firstErr != nil || nextIndex >= jobCount {
				mutex.Unlock()
				return
			}
			index := nextIndex
			nextIndex++
			mutex.Unlock()

			if err := job(index); err != nil {
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

func runIndexJobsWithSharedWorkers(jobCount int, workerGate chan struct{}, job func(int) error) error {
	if jobCount == 0 {
		return nil
	}
	if cap(workerGate) < 1 {
		return fmt.Errorf("shared workers must be at least one")
	}

	var firstErr error
	var mutex sync.Mutex
	var waitGroup sync.WaitGroup

	for index := 0; index < jobCount; index++ {
		waitGroup.Add(1)
		go func(index int) {
			defer waitGroup.Done()

			workerGate <- struct{}{}
			defer func() {
				<-workerGate
			}()

			mutex.Lock()
			shouldStop := firstErr != nil
			mutex.Unlock()
			if shouldStop {
				return
			}

			if err := job(index); err != nil {
				mutex.Lock()
				if firstErr == nil {
					firstErr = err
				}
				mutex.Unlock()
			}
		}(index)
	}

	waitGroup.Wait()
	return firstErr
}
