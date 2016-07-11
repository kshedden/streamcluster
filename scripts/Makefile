path = /nfs/kshedden/Teal_Furnholm/StreamCluster/

.PHONY: tests small check clean_tests clean

clean_tests:
	rm -f $(path)Tests/sim1*
	rm -f $(path)Tests/small1*

clean: clean_tests

small:
	go run streamclust.go -basepath=$(path) -kmer=3 -nclust=10 -skip=10 -normalize -stem=small1 -fasta=NCBI_MOCK_PROT_SMALL_SORT.fasta.gz

check:
	go run check.go -basepath=$(path) -stem=small1 -fasta=NCBI_MOCK_PROT_SMALL_SORT.fasta.gz

tests:
	go run streamclust.go -basepath=$(path)Tests -kmer=3 -nclust=2 -skip=1 -stem=test1 -fasta=test1.fasta.gz
	go run gen_simdat.go -basepath=$(path)Tests -outstem=sim1 -slen=1000 -nseq=2000 -p1=0.1 -p2=0.1
	go run streamclust.go -basepath=$(path)Tests -kmer=3 -nclust=2 -skip=500 -stem=sim1 -fasta=sim1.fasta.gz
	go run test.go -basepath=$(path)Tests
