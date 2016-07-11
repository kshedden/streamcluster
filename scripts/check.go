package main

import (
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq/linear"
	"github.com/kshedden/streamcluster"
)

var (
	basepath string

	outstem string

	fasta_in string

	readmax int

	clusters []int64

	clustermap map[int64][]int64

	//   Query letter
	//    -  A  C  G  T
	// -  0 -5 -5 -5 -5
	// A -5 10 -3 -1 -4
	// C -5 -3  9 -5  0
	// G -5 -1 -5  7 -3
	// T -5 -4  0 -3  8
	needle = align.NW{
		{0, -5, -5, -5, -5},
		{-5, 10, -3, -1, -4},
		{-5, -3, 9, -5, 0},
		{-5, -1, -5, 7, -3},
		{-5, -4, 0, -3, 8},
	}
)

func load_clusters() {

	fname := fmt.Sprintf("%s_clusters.json.gz", outstem)
	fname = path.Join(basepath, fname)
	fid, err := os.Open(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	rdr, err := gzip.NewReader(fid)
	if err != nil {
		panic(err)
	}
	defer rdr.Close()
	enc := json.NewDecoder(rdr)
	err = enc.Decode(&clusters)
	if err != nil {
		panic(err)
	}

	clustermap = make(map[int64][]int64)
	for k, v := range clusters {
		vec := clustermap[v]
		vec = append(vec, int64(k))
		clustermap[v] = vec
	}
}

// TODO: there must be a better way
func score(i1, j1 int, seql []string) float64 {

	t1 := &linear.Seq{}
	t1.Alpha = alphabet.DNAgapped
	t1.AppendLetters(alphabet.BytesToLetters([]byte(seql[i1]))...)
	t2 := &linear.Seq{}
	t2.Alpha = alphabet.DNAgapped
	t2.AppendLetters(alphabet.BytesToLetters([]byte(seql[j1]))...)

	aln, err := needle.Align(t1, t2)
	if err != nil {
		panic(err)
	}

	fa := align.Format(t1, t2, aln, '-')
	if err != nil {
		panic(err)
	}

	a1 := []byte(fmt.Sprintf("%s", fa[0]))
	a2 := []byte(fmt.Sprintf("%s", fa[1]))

	nc := 0
	for j := 0; j < len(a1); j++ {
		if a1[j] != a2[j] {
			nc++
		}
	}
	return float64(nc) / float64(len(a1))
}

func main() {

	// Read flags
	flag.StringVar(&basepath, "basepath", "", "all files are in this directory")
	flag.StringVar(&fasta_in, "fasta", "", "fasta input file")
	flag.StringVar(&outstem, "stem", "", "stem for output files")
	flag.IntVar(&readmax, "readmax", 10000, "maximum number of sequence records to read")
	flag.Parse()

	load_clusters()

	// Set up a fasta sequence reader
	fname := path.Join(basepath, fasta_in)
	fid, err := os.Open(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	rdr, err := gzip.NewReader(fid)
	if err != nil {
		panic(err)
	}
	defer rdr.Close()
	fr := streamcluster.NewFastaReader(rdr)

	// Read all the sequences
	var seql []string
	for idx := 0; idx < readmax; idx++ {
		_, seq := fr.Next()
		if seq == "" {
			break
		}
		seql = append(seql, seq)
	}

	fmt.Printf("Within-cluster scores:\n")
	nx := 0.0
	sx := 0.0
outer1:
	for _, v := range clustermap {
		for i1 := 0; i1 < len(v); i1++ {
			for j1 := 0; j1 < i1; j1++ {
				s := score(i1, j1, seql)
				sx += s
				fmt.Printf("%v\n", s)
				nx++
				if nx > 10 {
					break outer1
				}
			}
		}
	}
	fmt.Printf("  %v\n", sx/nx)

	fmt.Printf("Between-cluster scores:\n")
	nx = 0
	sx = 0
outer2:
	for k1, v1 := range clustermap {
		for k2, v2 := range clustermap {
			if k1 == k2 {
				continue
			}
			if len(v1) == 0 || len(v2) == 0 {
				continue
			}
			if v1[0] > int64(readmax) || v2[0] > int64(readmax) {
				continue
			}
			s := score(int(v1[0]), int(v2[0]), seql)
			sx += s
			fmt.Printf("%v\n", s)
			nx++
			if nx > 10 {
				break outer2
			}
		}
	}
	fmt.Printf("  %v\n", sx/nx)
}
