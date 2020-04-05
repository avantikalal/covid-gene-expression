#! /bin/bash

######################################################
# Extracting raw counts for DE analysis using DESeq2 #
######################################################

cd counts

# Obtaining the list of samples
ls -1 final | awk '{OFS="\t"; printf "%s\t%s/StringTie_transcripts_filtered.gtf\n",$1,$1}' > sample_list.txt
mv sample_list.txt final

cd final

# Extracting the counts for DE software using a script developed by StringTie
# Requires Python 2!!
python ../../scripts/prepDE.py -i sample_list.txt

echo "Count extraction done!"