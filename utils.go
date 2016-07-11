package streamcluster

import (
	"bufio"
	"io"
	"math"
	"strings"
)

type FastaReader struct {
	rdr     io.ReadCloser
	scanner *bufio.Scanner
	buf     []string
}

func NewFastaReader(rdr io.ReadCloser) *FastaReader {
	scanner := bufio.NewScanner(rdr)
	fr := &FastaReader{rdr, scanner, nil}
	fr.scan()
	return fr
}

func (fr *FastaReader) Next() (string, string) {

	fr.scan()

	if len(fr.buf) == 0 {
		return "", ""
	}

	if !strings.HasPrefix(fr.buf[0], ">") {
		panic("FastReader invalid state")
	}

	name := fr.buf[0][1:]
	fr.buf = fr.buf[1:]

	var sl []string
	for (len(fr.buf) > 0) && !strings.HasPrefix(fr.buf[0], ">") {
		sl = append(sl, fr.buf[0])
		fr.buf = fr.buf[1:]
	}

	return name, strings.Join(sl, "")
}

func (fr *FastaReader) scan() {
	// scan to next header or EOF
	for fr.scanner.Scan() {
		line := fr.scanner.Text()
		fr.buf = append(fr.buf, line)
		if strings.HasPrefix(line, ">") {
			return
		}
	}
	return
}

func GenHash(seq string, kmer_length int, normalize bool) []float64 {

	mx := uint64(math.Pow(4, float64(kmer_length)))
	fx := mx / uint64(4)

	count := make([]float64, mx)

	// Convert the sequence to integers
	iseq := make([]uint64, len(seq))
	for i, v := range seq {
		switch v {
		case 'A':
			iseq[i] = 0
		case 'T':
			iseq[i] = 1
		case 'G':
			iseq[i] = 2
		case 'C':
			iseq[i] = 3
		}
	}

	// Initialize the hash
	hash := uint64(0)
	b := uint64(1)
	for i := 0; i < kmer_length; i++ {
		hash += b * iseq[i]
		b *= 4
	}
	count[hash] += 1.0

	// Roll through the sequence
	for i := kmer_length; i < len(seq); i++ {
		hash /= 4
		hash += uint64(fx) * iseq[i]
		count[hash] += 1
	}

	if normalize {
		x := 0.0
		for _, y := range count {
			x += y
		}
		for i, _ := range count {
			count[i] /= x
		}
	}

	return count
}
