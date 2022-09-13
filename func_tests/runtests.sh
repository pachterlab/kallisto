#!/bin/bash

### Setup ###

kallisto="./src/kallisto"
bustools="bustools"
test_dir="./func_tests"

cmdexec() {
	cmd="$1"
	expect_fail="$2"
	printf "$cmd\n"
	if eval "$cmd 2> /dev/null 1> /dev/null"; then
		if [ "$expect_fail" = "1" ]; then
			printf "^[Failed - Returned normally but expected failure]\n"
			exit 1
		else
			printf "^[OK]\n"
		fi
	else
		if [ "$expect_fail" = "1" ]; then
			printf "^[OK - Exited non-zero as expected]\n"
		else
			printf "^[Failed]\n"
			exit 1
		fi
	fi
}

checkcmdoutput() {
	cmd="$1"
	correct_md5="$2"
	printf "$cmd\n"
	output_md5=$(eval "$cmd"|md5sum|awk '{ print $1 }')
	if [ "$output_md5" = "$correct_md5" ]; then
		printf "^[Output OK]\n"
	else
		printf "^[Output incorrect! Expected: "
		printf "$correct_md5"
		printf " Actual: "
		printf "$output_md5"
		printf "]\n"
		exit 1
	fi	
}

### FASTA files ###

echo ">t1
ACGTGATGAGTGAGTCAGT
>t2
ACGTGATGAGATGATGAGTCAGT
>t3
ACGTGATGTGAGTCAGT
>t4
ccccaaaaaa
>t5
ttttttgggg" > $test_dir/simple.fasta

echo ">t1
AAANNNTTTKK
>t2
NNNCCCCNCCC" > $test_dir/nonATCG.fasta

echo ">t1
AAAAAAAAAAAAAAAA
>t2
CCCAAAAAAAAAAAAA
>t3
CCCGGAAAAAAAAAAA
>t4
AAAAAAAAAAAAATTT
>t5
TAAAAAAAAAAAAAAT" > $test_dir/polyA.fasta

echo ">t1
AAATTTCCCGGGAAATTTCCCGGG
>t2
AAAtttCCCgggAAAtttCCCGGG
>t3
AAATTTCCCGGGAAATTTCCCGGG
>t4
AAATTTCCCGGGAAATTTCCCGGG
>t4
CCCCCCCCCCCCCCCCCCCCCCCC
>t5
CCCCCCCCCCCCCCCCCCCCCCCC
>t5
CCCCCCCCCCCCCCCCCCCCCCCC" > $test_dir/duplicates.fasta

### FASTQ files ###

cat $test_dir/simple.fasta $test_dir/nonATCG.fasta $test_dir/polyA.fasta $test_dir/duplicates.fasta | \
	awk '{if(NR%2){ print "\@" $0 } else { printf "%s\n+\n",$0; i=0; while (i++ < length($0)) printf "J"; printf "\n" } }' \
	| gzip > $test_dir/small.fastq.gz


cat /dev/null > $test_dir/medium.fastq.gz
for i in {1..700}; do cat $test_dir/small.fastq.gz >> $test_dir/medium.fastq.gz; done

cat /dev/null > $test_dir/large.fastq.gz
for i in {1..100}; do cat $test_dir/medium.fastq.gz >> $test_dir/large.fastq.gz; done

rm $test_dir/medium.fastq.gz

# Paired-end reads for basic7.idx index

echo "@t1
ACGTGATG
+
JJJJJJJJ
@t2
ACGTGATG
+
JJJJJJJJ
@t3
CCCCCCCC
+
JJJJJJJJ
@t4
ccccaaaaaa
+
JJJJJJJJJJ
@t5
ttttttgggg
+
JJJJJJJJJJ
" | gzip > $test_dir/simple_pair1.fastq.gz

echo "@t1
GAGTCAGT
+
JJJJJJJJ
@t2
TGAGTCAG
+
JJJJJJJJ
@t3
ACGTGATG
+
JJJJJJJJ
@t4
ccccaaaaaa
+
JJJJJJJJJJ
@t5
ccccaaaaaa
+
JJJJJJJJJJ
" | gzip > $test_dir/simple_pair2.fastq.gz

# Barcode+UMI FastQ file for 10Xv3

echo "@t1
GCCCCCCCCCCCCCCGTAAAAAAAAAATGG
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@t2
GCCCCCCCCCCCCCCGTAAAAAAAAAATGG
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@t3
GCCCCCCCCCCCCCCGTAAAAttAAAAtGG
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@t4
GCCCCCCggCCCCCCGtAAAAttAAAAtGG
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@t5
GCCCCCCttCCCCCCGtAAAAggAAAAtGG
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
" | gzip > $test_dir/10xv3.fastq.gz


### TEST - kallisto index ###

echo "[INDEX]"

# Test k=7

cmdexec "$kallisto index -i $test_dir/basic7.idx -k 7 $test_dir/simple.fasta"

# Test k=9

cmdexec "$kallisto index -i $test_dir/basic9.idx -k 9 $test_dir/simple.fasta"

# Test non-ATCG bases

cmdexec "$kallisto index -i $test_dir/nonATCG.idx -k 5 $test_dir/nonATCG.fasta"

# Test polyA truncation

cmdexec "$kallisto index -i $test_dir/polyA.idx -k 5 $test_dir/polyA.fasta"

# Test k = even number (should fail)

cmdexec "$kallisto index -i $test_dir/basic_fail.idx -k 10 $test_dir/simple.fasta" 1

# Test duplicate transcript names (should fail)

cmdexec "$kallisto index -i $test_dir/duplicates_fail.idx -k 11 $test_dir/duplicates.fasta" 1

# Test duplicate transcript names with --make-unique

cmdexec "$kallisto index -i $test_dir/duplicates.idx -k 11 --make-unique $test_dir/duplicates.fasta"


### TEST - kallisto quant ###

# Test single-end small.fastq.gz using various indices

cmdexec "$kallisto quant -o $test_dir/quantbasic -i $test_dir/basic7.idx --single -l 5 -s 2 $test_dir/small.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasic/abundance.tsv" f1c8927dc29b943758242902d5c45a86

cmdexec "$kallisto quant -o $test_dir/quantnonATCG -i $test_dir/nonATCG.idx --single -l 5 -s 2 $test_dir/small.fastq.gz"
checkcmdoutput "cat $test_dir/quantnonATCG/abundance.tsv" 9272185ff011f9a840c29c5c40a2d390

cmdexec "$kallisto quant -o $test_dir/quantpolyA -i $test_dir/polyA.idx --single -l 5 -s 2 $test_dir/small.fastq.gz"
checkcmdoutput "cat $test_dir/quantpolyA/abundance.tsv" 745539c18d4ff07bf6de0938563f2362

cmdexec "$kallisto quant -o $test_dir/quantduplicates -i $test_dir/duplicates.idx --single -l 5 -s 2 $test_dir/small.fastq.gz"
checkcmdoutput "cat $test_dir/quantduplicates/abundance.tsv" 49a62b913a155598dd6894a2f0eb0905

# Test single-end small.fastq.gz strand-specificity

cmdexec "$kallisto quant -o $test_dir/quantbasicfr -i $test_dir/basic7.idx --single --fr-stranded -l 5 -s 2 $test_dir/small.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasicfr/abundance.tsv" ce2ed5a3a1bab582fcb62dc02f4d9323

cmdexec "$kallisto quant -o $test_dir/quantbasicrf -i $test_dir/basic7.idx --single --rf-stranded -l 5 -s 2 $test_dir/small.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasicrf/abundance.tsv" 017e8ba77d7e7b39a60bb7c047e620dc

# Test paired-end small.fastq.gz

cmdexec "$kallisto quant -o $test_dir/quantbasicpaired -i $test_dir/basic7.idx $test_dir/simple_pair1.fastq.gz $test_dir/simple_pair2.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasicpaired/abundance.tsv" 1cf9cd508c04eff7cbda6e74cc1e44ea

# Test paired-end small.fastq.gz with multiple pairs of files and with strand-specificity

cmdexec "$kallisto quant -o $test_dir/quantbasicpairedfr -i $test_dir/basic7.idx --fr-stranded $test_dir/simple_pair1.fastq.gz $test_dir/simple_pair2.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasicpairedfr/abundance.tsv" bc6a75291c2eb661a57b55de534dd8f0

cmdexec "$kallisto quant -o $test_dir/quantbasicpairedrf -i $test_dir/basic7.idx --rf-stranded $test_dir/simple_pair1.fastq.gz $test_dir/simple_pair2.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasicpairedrf/abundance.tsv" dd1ed6ed8fb393616817a3dd506f821f

cmdexec "$kallisto quant -o $test_dir/quantbasicpairedmultfr -i $test_dir/basic7.idx --fr-stranded $test_dir/simple_pair1.fastq.gz $test_dir/simple_pair2.fastq.gz $test_dir/simple_pair2.fastq.gz $test_dir/simple_pair1.fastq.gz"
checkcmdoutput "cat $test_dir/quantbasicpairedmultfr/abundance.tsv" f810d28aed3969514f1d21120a6fc825

# Test multiple large fastq files with more threads

cmdexec "$kallisto quant -o $test_dir/quantlarge -t 12 -i $test_dir/basic7.idx --single -l 5 -s 2 $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz"
checkcmdoutput "cat $test_dir/quantlarge/abundance.tsv" 024e56c50c774aa09e00a441512415aa


### TEST - kallisto bus ###

if ! command -v bustools &> /dev/null
then
    echo "Error: bustools could not be found"
    exit 1
fi

# Try custom "magic string" for -x with multiple large files and many threads

cmdexec "$kallisto bus -o $test_dir/buslarge -t 12 -i $test_dir/basic7.idx -x 0,0,2:0,2,6:1,0,0 $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz $test_dir/large.fastq.gz"
cmdexec "$bustools sort -o $test_dir/buslarge/output.s.bus -t 12 $test_dir/buslarge/output.bus"
checkcmdoutput "bustools text -p $test_dir/buslarge/output.s.bus|cut -f1,2,4" 2b4e7120a9ee3d419c1de5e1689d1634

# Try a simple 10XV3 (with unstranded pseudoalignment)

cmdexec "$kallisto bus -o $test_dir/bus10xv3 -t 1 -i $test_dir/basic7.idx -x 10XV3 --unstranded $test_dir/10xv3.fastq.gz $test_dir/small.fastq.gz"
cmdexec "$bustools sort -o $test_dir/bus10xv3/output.s.bus -t 12 $test_dir/bus10xv3/output.bus"
checkcmdoutput "$bustools text -p $test_dir/bus10xv3/output.s.bus" f8f76f64d3ef351057f3bae016e2fa8a


# Try processing demultiplexed bulk RNA-seq with strand-specificity with EM and kallisto quant-tcc (and compare with quant) 

echo "t1	g1
t2	g1
t3	g2
t4	g2
t5	g3" > $test_dir/t2g_test.txt # Make a transcript-to-gene mapping file
cmdexec "$kallisto bus --fr-stranded -o $test_dir/busbulklarge -t 12 -i $test_dir/basic7.idx $test_dir/large.fastq.gz"
cmdexec "$bustools sort -o $test_dir/busbulklarge/output.s.bus -t 12 $test_dir/busbulklarge/output.bus"
cmdexec "$bustools count --cm -m -o $test_dir/busbulklarge/counts_tcc/ -g $test_dir/t2g_test.txt -t $test_dir/busbulklarge/transcripts.txt -e $test_dir/busbulklarge/matrix.ec $test_dir/busbulklarge/output.bus"
cmdexec "$kallisto quant-tcc -o $test_dir/quant_tcc_test1/ -l 5 -s 2 -i $test_dir/busbulklarge/index.saved -g $test_dir/t2g_test.txt -e $test_dir/busbulklarge/counts_tcc/output.ec.txt $test_dir/busbulklarge/counts_tcc/output.mtx"
checkcmdoutput "cat $test_dir/quant_tcc_test1/matrix.abundance.gene.tpm.mtx" ed7c7fa08e283cc6a1b136cc0e1ab039



