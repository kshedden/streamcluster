package main

// Simulate idealized two-group data for testing

import (
	"compress/gzip"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"path"
)

var (
	//  All files are located in this directory
	basepath string

	// Stem of all output files
	outstem string

	// Number of sequences to simulate
	nseq int

	// Length of each simulated sequence
	slen int

	// Mutation probability for first cluster
	pmiss1 float64

	// Mutation probability for second cluster
	pmiss2 float64
)

// Generate a random letter, uniformly from ATGC
func gen_letter() byte {
	x := rand.Int63n(4)
	switch x {
	case 0:
		return 'A'
	case 1:
		return 'T'
	case 2:
		return 'G'
	case 3:
		return 'C'
	}
	return 'A' // Cannot reach here
}

// Generate a random sequence of ATGC letters with the given length
func gen_seq(slen int) []byte {
	iseq := make([]byte, slen)
	for i := 0; i < slen; i++ {
		iseq[i] = gen_letter()
	}
	return iseq
}

func main() {

	// Read flags
	flag.StringVar(&basepath, "basepath", "", "all files are in this directory")
	flag.StringVar(&outstem, "outstem", "", "stem for output file names")
	flag.IntVar(&nseq, "nseq", 20, "number of sequences to generate")
	flag.IntVar(&slen, "slen", 1000, "length of each sequence")
	flag.Float64Var(&pmiss1, "p1", 0, "")
	flag.Float64Var(&pmiss2, "p2", 0, "length of each sequence")
	flag.Parse()

	// Cluster centers
	seq1 := gen_seq(slen)
	seq2 := gen_seq(slen)

	fname := fmt.Sprintf("%s.fasta.gz", outstem)
	fname = path.Join(basepath, fname)
	fid, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	wtr := gzip.NewWriter(fid)
	defer wtr.Close()

	for k := 0; k < nseq; k++ {
		wtr.Write([]byte(fmt.Sprintf(">S%03d\n", k)))

		var pmiss float64
		var seq []byte
		if k < nseq/2 {
			pmiss = pmiss1
			seq = seq1
		} else {
			pmiss = pmiss2
			seq = seq2
		}

		iseq1 := make([]byte, slen)
		for j := 0; j < slen; j++ {
			if rand.Float64() < pmiss {
				iseq1[j] = gen_letter()
			} else {
				iseq1[j] = seq[j]
			}
		}

		for len(iseq1) > 0 {
			jj := 70 // Line length
			if len(iseq1) < jj {
				jj = len(iseq1)
			}
			wtr.Write(iseq1[0:jj])
			wtr.Write([]byte("\n"))
			iseq1 = iseq1[jj:]
		}
	}
}
