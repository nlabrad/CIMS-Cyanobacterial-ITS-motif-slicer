#ITS seq processing with BASH 

#basic grep to find seqs and make a variable
#myVar=$(grep -o 'GACCT[ATCGR]*AGGTC' GeitTest1.txt)
#echo ${myVar} > FILENAME.txt
#if you want to save into one file with all seqs
#for file in *.txt; do; boxBseq=$(grep -o 'AGCAC[ATCGR]*ATGCT' $file); echo $boxBseq > ./boxB/boxBs; done;


#1 - removes lines starting with > and saves each sepearately then removes spaces from files and saves to new file 
for file in *.fasta; do; grep -v '>' $file | perl -pe 's/\s+//g' > $file.txt; done;

#2 - finds D1D1 and puts into D1 folder
mkdir D1D1seqs
for file in *.txt; do; D1D1seq=$(grep -o 'GACCT[ATCGR]*AGGTCA' $file); echo $D1D1seq > ./D1D1seqs/D1D1$file; done;
#** this is going to the second aggtc in the seq not the first - had to add A 

#3 - finds boxB
mkdir boxB
for file in *.txt; do; boxBseq=$(grep -o 'AGCAC[ATCGR]*ATGCT' $file); echo $boxBseq > ./boxB/boxB$file; done;

#4 - finds tRNA genes
#tRNA-Ile always  begins with GGGCTATTA and ends with GGCCCA (74 bp long)
#tRNA-Ala always begins with GGGG and ends with CCTCCA (73 bp long)
for file in *.fasta; do; tRNA_Ile=$(grep -o 'GGGCTATTA[ATCGR]*GGCCCA' $file); echo $tRNA_Ile > tRNA_Iles.txt; done;
for file in *.fasta; do; tRNA_Ala=$(grep -o 'GGGG[ATCGR]*CCTCCA' $file); echo $tRNA_Ala > tRNA_Ala.txt; done;

#v3
start - GUAGG
end - CAGAC



How to make sure motif is the right size 

for file in D1*; do;
    if [$(wc -m $file) == 0];
    then 
        rm $file
    else 
        echo "This region exists"
    fi;
done;




notes: 
to remove the title line and all spaces between lines  -- not sure if works yet
grep -v '>' FILE.txt | perl -pe 's/\s+//g' 

