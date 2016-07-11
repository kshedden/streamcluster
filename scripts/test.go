package main

import (
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path"
)

var (
	basepath string
)

func test(stem string, expclust []int64) {

	fname := fmt.Sprintf("%s_clusters.json.gz", stem)
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
	dec := json.NewDecoder(rdr)
	var clust []int64
	err = dec.Decode(&clust)
	if err != nil {
		panic(err)
	}

	for i, x := range clust {
		if x != expclust[i] {
			panic("Fail")
		}
	}
}

func main() {

	// Read flags
	flag.StringVar(&basepath, "basepath", "", "all files are in this directory")
	flag.Parse()

	test1_json := []byte("[0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1]")
	var expclust []int64
	json.Unmarshal(test1_json, &expclust)
	test("test1", expclust)
}
