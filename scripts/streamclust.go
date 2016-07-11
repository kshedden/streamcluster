package main

// Cluster a collection of sequences in a fasta sequence file using
// streaming k-means.
//
// Rough outline of the algorithm:
//
// 1. Randomly select m sequences to serve as the initial cluster centers.
// 2. Cycle through the sequences as they are read from the file, assign
//    each sequence to its closest centroid, then update the mean of that
//    centroid.
//
// Only one pass is made through the data.

import (
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"path"
	"sync"

	"github.com/kshedden/streamcluster"
)

var (
	// All files are in this directory
	basepath string

	// The hashes are based on k-mers of this length
	kmer_length int

	// Number of clusters
	nclust int

	// Randomly skip from skip to 2*skip sequences to initialize
	// cluster centers
	skip int

	// Stem for all output files
	outstem string

	// Sequence data file
	fasta_in string

	// Cluster centers
	centers [][]float64

	// Cluster sizes
	csize []float64

	// Cluster assignments
	clusters []int

	// If true, kmer distributions are normalized
	normalize bool

	// Need to protect the centroids while updating
	mutex sync.Mutex

	wg sync.WaitGroup

	xchan chan *xrec_t

	// Some channels used to manage the concurrency.
	done1    chan bool
	done2    chan bool
	throttle chan bool
)

// Randomly select sequences to serve as the initial cluster centers
func init_seeds(nseed int) [][]float64 {

	fmt.Printf("Getting seeds...")

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

	var hashes [][]float64
	var seq string

	idx := -1
	for j := 0; j < nseed; j++ {

		if j%100 == 0 {
			fmt.Printf(".")
		}

		// Skip a random number of sequences
		jmp := rand.Int63n(int64(skip)) + int64(skip)
		for k := 0; k < int(jmp); k++ {
			_, seq = fr.Next()
			idx++
		}

		hash := streamcluster.GenHash(seq, kmer_length, normalize)

		// Convert to float
		fhash := make([]float64, len(hash))
		for k, x := range hash {
			fhash[k] = float64(x)
		}

		hashes = append(hashes, fhash)
	}

	fmt.Printf("done\n")

	return hashes
}

// Returns the distance from the given hash to cluster center i.
func distance(hash []float64, i int) float64 {

	var d float64

	for j, x := range hash {
		mutex.Lock()
		u := float64(x) - centers[i][j]
		mutex.Unlock()
		d += u * u
	}

	return d
}

// The argmin of the given float slice.
func minloc(x []float64) int {

	i := 0
	v := x[0]

	for j := 1; j < len(x); j++ {
		if x[j] < v {
			v = x[j]
			i = j
		}
	}

	return i
}

type xrec_t struct {
	idx  int
	iloc int
	hash []float64
}

// The core clustering calculations, can be run concurrently.
func clustcore(seq string, idx int) {

	defer func() { <-throttle }()
	defer wg.Done()

	distances := make([]float64, nclust)
	hash := streamcluster.GenHash(seq, kmer_length, normalize)

	// Get the cluster assignment
	for i := 0; i < nclust; i++ {
		distances[i] = distance(hash, i)
	}
	iloc := minloc(distances)

	xchan <- &xrec_t{idx, iloc, hash}
}

func stream() {
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

	fmt.Printf("Clustering...")

	// Process the sequences
	go func() {
		for idx := 0; ; idx++ {
			_, seq := fr.Next()
			if seq == "" {
				done1 <- true
				return
			}

			wg.Add(1)
			throttle <- true
			go clustcore(seq, idx)
		}
	}()

	// Retrieve the results
	go func() {
		nd := 0
		for r := range xchan {

			// Update the progress meter
			nd++
			if nd%100 == 0 {
				fmt.Printf(".")
				if nd%1000 == 0 {
					fmt.Printf("[%d]", nd)
				}
			}

			// Update the cluster center
			mutex.Lock()
			c := centers[r.iloc]
			for j := 0; j < len(centers[r.iloc]); j++ {
				c[j] *= csize[r.iloc]
				c[j] += r.hash[j]
				c[j] /= csize[r.iloc] + 1
			}
			mutex.Unlock()
			csize[r.iloc]++

			// Update the cluster membership
			if len(clusters) < r.idx+1 {
				// Add space if needed
				clusters = append(clusters, make([]int, r.idx+1-len(clusters))...)
			}
			clusters[r.idx] = r.iloc
		}
		done2 <- true
	}()

	<-done1   // Done reading file
	wg.Wait() // Done processing sequences
	close(xchan)
	<-done2 // Done updating centroids and cluster assignments

	fmt.Printf(" done\n")
}

func write_results() {

	fname := fmt.Sprintf("%s_clusters.json.gz", outstem)
	fname = path.Join(basepath, fname)
	fid, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	wtr := gzip.NewWriter(fid)
	defer wtr.Close()
	enc := json.NewEncoder(wtr)
	err = enc.Encode(clusters)
	if err != nil {
		panic(err)
	}

	fname = fmt.Sprintf("%s_centers.json.gz", outstem)
	fname = path.Join(basepath, fname)
	fid, err = os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	wtr = gzip.NewWriter(fid)
	defer wtr.Close()
	enc = json.NewEncoder(wtr)
	err = enc.Encode(centers)
	if err != nil {
		panic(err)
	}
}

func main() {

	// Read flags
	flag.StringVar(&basepath, "basepath", "", "all files are in this directory")
	flag.StringVar(&fasta_in, "fasta", "", "fasta input file")
	flag.IntVar(&kmer_length, "kmer", 3, "k-mer length")
	flag.IntVar(&nclust, "nclust", 100, "number of clusters")
	flag.IntVar(&skip, "skip", 100, "skip sequences during initialization")
	flag.BoolVar(&normalize, "normalize", false, "normalize kmer profiles")
	flag.StringVar(&outstem, "stem", "", "stem for output files")
	flag.Parse()

	centers = init_seeds(nclust)
	csize = make([]float64, nclust)

	xchan = make(chan *xrec_t)

	// The size of this channel is the maximum number of sequences
	// that are concurrently being processed.
	throttle = make(chan bool, 500)

	done1 = make(chan bool, 1)
	done2 = make(chan bool, 1)

	stream()
	write_results()
}
