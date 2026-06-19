package app

import (
	workerpool "atsp_aco_msa/modules/experiments/workers"
	"fmt"
	"sort"
)

func instanceJobsByDescendingDimension(atspsData []AtspData, label string, run func(AtspData) error) []workerpool.Job {
	orderedAtspData := append([]AtspData{}, atspsData...)
	sort.SliceStable(orderedAtspData, func(i, j int) bool {
		leftDimension := len(orderedAtspData[i].Matrix)
		rightDimension := len(orderedAtspData[j].Matrix)
		if leftDimension != rightDimension {
			return leftDimension > rightDimension
		}

		return orderedAtspData[i].Name < orderedAtspData[j].Name
	})

	jobs := make([]workerpool.Job, 0, len(orderedAtspData))
	for _, atspData := range orderedAtspData {
		atspData := atspData
		jobs = append(jobs, workerpool.Job{
			Label: fmt.Sprintf("%s/%s", label, atspData.Name),
			Run: func() error {
				return run(atspData)
			},
		})
	}

	return jobs
}
