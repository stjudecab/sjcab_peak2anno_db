
Blacklists from various sources might be slightly different. Here on 2023-04-11 we merge all previously scattered sources of blacklists and encourage everyone to use this merged collection.
The main sources are:
1. Blacklists from Boyle lab: https://github.com/Boyle-Lab/Blacklist/tree/master/lists
2. Blacklists from Encode: https://www.encodeproject.org/annotations/ENCSR636HFF/
3. Blacklists from Anshul Kundaje: https://personal.broadinstitute.org/anshul/projects/encode/rawdata/blacklists/hg19-blacklist-README.pdf
4. for sacCer3 the blacklist was done in-house, based on the Input samples from CAB-1112 project, task CAB-1206; as described here in details: https://www.evernote.com/l/Aju0gjqFnKFOBK6GyZ085aBVqv8o0PqYlN4

Commands to merge blacklists that we had prior April 2023:

```
for genome in ce10 ce11 dm3 dm6 hg18 hg19 hg38 mm10 mm39 mm9 sacCer3
do
    for bkFile in $(tree -i -f ../../.. | grep "$genome-blacklist")
    do
        awk '{print $1 "\t" $2 "\t" $3}' $bkFile >> tmp.bed
    done
    bedSort tmp.bed tmp.bed
    bedtools merge -i tmp.bed > $genome-blacklist.bed
    rm tmp.bed
done
wc -l *.bed > stats.txt
```

Stats are:

```
   100 ce10-blacklist.bed
    97 ce11-blacklist.bed
   271 dm3-blacklist.bed
   182 dm6-blacklist.bed
  1411 hg18-blacklist.bed
  1871 hg19-blacklist.bed
  2375 hg38-blacklist.bed
  4874 mm10-blacklist.bed
  3004 mm39-blacklist.bed
  3038 mm9-blacklist.bed
    54 sacCer3-blacklist.bed
```

Note, that if the collection of blacklists would be updated in the future, the files here should be the first ones to be updated.

