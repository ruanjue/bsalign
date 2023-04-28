# bsalign

bsalign is a library/tool for adaptive banding striped 8/2-bit-scoring global/extend/overlap DNA sequence pairwise/multiple alignment

# Installation
```sh
git clone https://github.com/ruanjue/bsalign.git
cd bsalign
make
```

# run bsalign

```sh
bsalign
```
```txt
commands:
 align       Pairwise alignment implemented by 8-bit encoded Banded Striped SIMD
 edit        Pairwise alignment using edit distance implemented by 2-bit encoded banded Striped algorithm
 poa         Multiple alignment implemented by 8-bit encoded Banded Striped SIMD Partial Order Alignment
 cat         Concatenate pieces of seqs into one seq by overlaping
```

## Example
```sh
cd example
sh run.sh
```

## Result Example
```txt
29.1	75	+	0	75	29.2	76	+	0	76	128	0.934	71	4	0	1
TGTTACTTTTCTTCCCTGCTGTATAAACCC-CAGTTTTAGTCAGTCAGGGAGATGGATTTGAGACTGAGCTCCCAT
||||||*|||||||||||||**||||||||-*||||||||||||||||||||||||||||||||||||||||||||
TGTTACATTTCTTCCCTGCTACATAAACCCTTAGTTTTAGTCAGTCAGGGAGATGGATTTGAGACTGAGCTCCCAT
```

## Result Format
Each result is 4 lines.
Line 1 :col1-RefName; col2-RefLength; col3-RefStrand; col4-RefStart; col5-RefEnd; col6-QueryName; col7-QueryLength; col8-QueryStrand; col9-QueryStart; col10-QueryEnd; col11-AlignmentScore; col12-Identity; col13-NumberOfMatch; col14-NumberOfMismatch; col15-NumberOfDeletion; col16-NumberOfInsertion
Line 2 :Reference Sequence
Line 3 :'|', '*' and '-' mean match, mismatch and indel, respectively.
Line 4 :Query Sequence

# use bsalign library

copy bsalign directory into your code `cp -r /path/to/bsalign .`

## Pairwise Alignment Example
```txt
#include "bsalign/bsalign.h"

int verbose = 0; // be quiet in alignment
b1i mtx[16]; // score matrix, 4 * 4
banded_striped_epi8_seqalign_set_score_matrix(mtx, sc_mat=2, sc_mis=-6); // init score matrix
b1v *memp = adv_init_b1v(1024, 0, WORDSIZE, 0); // it needs a WORDSIZE(16 bytes)-aligned memory block to perform SIMD alignment
u4v *cigars = init_u4v(32); // use to store alignment cigars (SAM-like), or NULL if useless
int bandwidth = 128; // Or 0 if disable banded alignment
// perform pairwise global alignment (8-bits)
seqalign_result_t rs = banded_striped_epi8_seqalign_pairwise((u1i*)qseq, qlen, (u1i*)tseq, tlen, memp, cigars, SEQALIGN_MODE_GLOBAL, bandwidth, mtx, sc_gapo=-3, sc_gape=-2, 0, 0, verbose);
// perform pairwise global edit (2-bits), using edit-distance in alignment, much faster than 8-bits alignment
seqalign_result_t rs = striped_seqedit_pairwise((u1i*)qseq, qlen, (u1i*)tseq, tlen, SEQALIGN_MODE_GLOBAL, bandwidth, memp, cigars, verbose);
// perform pairwise kmer-guided edit (2-bits), it is better for two strange reads, because it infers the outline of alignment by kmer-matching-synteny
seqalign_result_t rs = kmer_striped_seqedit_pairwise(ksize=13, (u1i*)qseq, qlen, (u1i*)tseq, tlen, memp, cigars, verbose);
// print alignment information
fprintf(stdout, "QRY\t%d\t%d\tREF\t%d\t%d\tmat=%d\tmis=%d\tins=%d\tdel=%d\n", rs.qb, rs.qe, rs.tb, rs.te, rs.mat, rs.mis, rs.ins, rs.del);
char *alnstr[3];
alnstr[0] = malloc(rs.aln + 1);
alnstr[1] = malloc(rs.aln + 1);
alnstr[2] = malloc(rs.aln + 1);
seqalign_cigar2alnstr(qseq, tseq, &rs, cigars, alnstr, rs.aln);
// print alignment string
fprintf(stdout, "%s\n%s\n%s\n", alnstr[0], alnstr[2], alnstr[1]);
free(alnstr[0]); free(alnstr[1]); free(alnstr[2]);
free_u4v(cigars);
free_b1v(memp);
```

## Multiple Alignment Example
```txt
#include "bsalign/bspoa.h"

BSPOAPar par = DEFAULT_BSPOA_PAR; // change par.xxx if you want
BSPOA *g = init_bspoa(par);
beg_bspoa(g); // prepare to accept reads
for(...) push_bspoa(g, (char*)rdseq, (int)rdlen); // push reads one by one
end_bspoa(g); // MSA generated
tidy_msa_bspoa(g); // polish MSA to call more SNVs
call_snvs_bspoa(g); // call SNVs on the polished MSA
// print MSA, linewidth=0 to output each read in a single line
// colorful=1 to output friendly terminal characters, pipe to 'less -S -R' if no color in your screen
print_msa_bspoa(g, "<MSA_ID>", 0, 0, linewidth=100, colorful=1, stdout);
print_snvs_bspoa(g, "<MSA_ID>", stdout);
// Or write binary MSA (no SNVs) to save disk space
dump_binary_msa_bspoa(g, "Welcome to AGIS", 15, file);
// Load a binary MSA instead of beg/push/end_bspoa, Note: invoke call_snvs_bspoa if you want SNVs
String *metainfo = init_string(32);
load_binary_msa_bspoa(g, file, metainfo);
free_bspoa(g);
```
## Operating system
The Linux platform is supported by BSAlign. It run and tested on Ubuntu 20.04.2 LTS and CentOS Linux 7.

# Contact
Jue Ruan <ruanjue@gmail.com> <br>
Jue Ruan <ruanjue@caas.cn>
