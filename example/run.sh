FA=real.ont.b10M.txt
# parameter -M 2 -X 2 -O 4 -E 2 -Q 0 -P 0

echo Running NoBand.bsalign
../bsalign align -M 2 -X 2 -O 4 -E 2 -Q 0 -P 0 $FA > $FA.NoBand.bsalign 
echo Running Band64.bsalign
../bsalign align -W 64 -M 2 -X 2 -O 4 -E 2 -Q 0 -P 0 -m overlap  $FA >$FA.Band64.bsalign
echo Running Edit.bsalign
../bsalign edit -W 0 $FA > $FA.Edit0.bsalign
