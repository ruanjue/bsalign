# bsalign

bsalign is a library/tool for adaptive banding striped 8-bit-scoring global/extend/overlap DNA sequence alignment

# Installation
```sh
git clone https://github.com/ruanjue/bsalign.git
cd bsalign
make
```

# run bsalign

```sh
bslaign -h
```
### options
```txt
 -h          Show this document
 -m <string> align mode: global, extend, overlap [overlap]
 -W <int>    Bandwidth, 0: full length of query [0]
             set it to a small value if the begins of two reads are nearly aligned
 -M <int>    Score for match, [2]
 -X <int>    Penalty for mismatch, [6]
 -O <int>    Penalty for gap open, [3]
 -E <int>    Penalty for gap extension, [2]
 -Q <int>    Penalty for gap2 open, [0]
 -P <int>    Penalty for gap2 extension, [0]
 -R <int>    repeat times (for benchmarking) [1]
 -v          Verbose
 Gap-weighting functions
 If P < E and Q + P > O + E, 2-piecewise affine gap cost
  e.g. -M 2 -X 6 -O 3 -E 2 -Q 18 -P 1
 Else if O > 0, 1-piecewise affine gap cost
  e.g. -M 2 -X 6 -O 3 -E 2 -Q 0 -P 0
 Else, linear gap cost
  e.g. -M 2 -X 6 -O 0 -E 5 -Q 0 -P 0
```

# Contact
Jue Ruan <ruanjue@gmail.com> <br>
Jue Ruan <ruanjue@caas.cn>
