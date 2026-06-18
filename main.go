package main

import (
	"atsp_aco_msa/modules/app"
	"os"
)

func main() {
	app.Run(os.Args[1:])
}
