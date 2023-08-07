#!/bin/bash
if [ -z "$1" ]; then
  echo 'No directory supplied' >&2
  exit 1
fi

project=$1
now=$(date)
nowformat=$(date +"%y-%m-%d")

echo "-- create output folder $project.$nowformat; add first files"
mkdir -p $project.$nowformat

echo "-- moving in data dir"
cd $project

echo "-- loading config"
source config.txt

echo "-- creating directories"
mkdir -p logs
mkdir -p raw
mkdir -p tmp

if [ $classificationOnly -ne 1 ]
  then
if [ $skip_preprocessing -ne 1 ]
  then

  echo "-- decompressing raw files"
  #extracting files
  find . -name '*.gz' -print0 | xargs -0 -I {} -P $threads gunzip {}

  # looping through all files for merging
  echo "-- starting merging and filter "
  for f in *_R1_*.fastq; do

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(echo $f | sed "s/_L001.*//g")
  	total=$(grep "^@M" $f | wc -l)

    echo " "
    echo "===================================="
    echo "Processing sample $s"
    echo "(F: $f R: $r)"
    echo "===================================="
    $vsearch --fastq_mergepairs  $f \
          --reverse $r \
          --fastq_minovlen $merge_minovlen \
          --fastq_maxdiffs $merge_maxdiffs \
          --fastqout $s.merge.fq \
          --fastq_eeout \
          --relabel R1+2-${s}_ \
          --threads $threads  2> logs/vsearch.m.$s.log

    var="$(grep "Merged" logs/vsearch.m.$s.log)"
    echo "$s : $var" | tee -a logs/_merging.log

    $vsearch --fastq_filter $s.merge.fq \
      --fastq_stripleft $stripleft \
      --fastq_stripright $stripright \
      --fastq_maxee $fastq_maxee  \
      --fastq_minlen $fastq_minlen  \
      --fastq_maxlen $fastq_maxlen  \
      --fastq_maxns $fastq_maxns \
      --fastaout $s.mergefiltered.fa  \
      --fasta_width 0 \
      --threads $threads 2> logs/vsearch.mf.$s.log

    var="$(grep "sequences kept" logs/vsearch.mf.$s.log)"
    echo "$s : $var" | tee -a logs/_filter.log

    $vsearch --fastq_truncee $fastq_truncee \
          --fastq_filter $f \
          --fastq_minlen $fastq_minlen \
          --fastaout $s.trunc.fa \
          --relabel R1-${s}_ \
          --threads $threads 2> logs/vsearch.tf.$s.log

    var="$(grep "sequences kept" logs/vsearch.tf.$s.log)"
    echo "$s : $var" | tee -a logs/_truncfilter.log

  done

  cat *mergefiltered.fa > all.merge.fasta
  cat *trunc.fa > all.trunc.fasta

  if [ $use_fw_only -eq 1 ]
    then
      cp all.merge.fasta all.merge.fasta.bak
      cp all.trunc.fasta all.merge.fasta
      echo " "
      echo "!!! INFO: using forward reads only"
      echo " "
  fi #end useFWonly

    # cleanup
    echo "-- cleanup"
    mv *.fastq raw/
    mv *.fq tmp/
    mv *.fa tmp/

fi #end skippp

  echo " "
  echo "===================================="
  echo "ASV generation and mapping"
  echo "===================================="

  echo "-- derep"

  $vsearch --derep_fulllength all.merge.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.merge.derep.uc \
    --output all.merge.derep.fa \
    --threads $threads 2> logs/_derep.log

  echo "-- denoise"
  $vsearch --cluster_unoise  all.merge.derep.fa \
    --sizein --sizeout \
    --relabel ASV \
    --centroids zotus.merge_chim.fa \
    --threads $threads 2> logs/_unoise.log

  grep "Clusters" logs/_unoise.log
  grep "Singletons" logs/_unoise.log

  echo "-- sorting"
  $vsearch --sortbysize zotus.merge_chim.fa \
    --output zotus.merge_sorted.fa \
    --threads $threads 2>  logs/_sort.log

  echo "-- chimera removal"
  $vsearch --uchime3_denovo zotus.merge_sorted.fa \
    --abskew 16 \
    --nonchimeras zotus.merge.fa \
    --threads $threads 2>  logs/_uchime.log

  grep "Found" logs/_uchime.log

  ### create community table
  echo "-- add barcodelabel"
  cat all.merge.fasta |  sed "s/^>R1+2-\(.*\)\_\([0-9]*\);/>R1+2-\1_\2;barcodelabel=\1;/g" |  sed "s/^>R1-\([a-zA-Z0-9-]*\)\_\([0-9]*\)/>R1-\1_\2;barcodelabel=\1;/g" > all.merge.bc.fasta

  echo "-- map data against ASVs"
  $vsearch --usearch_global all.merge.bc.fasta --db zotus.merge.fa --strand plus --id 0.97 --uc map.merge.uc --otutabout asv.tab-csv --sizeout --threads $threads 2> logs/_mapping.log

  grep "Matching" logs/_mapping.log

  #python2.7 $p/uc2otutab.py map.trunc.uc > zotu_table.trunc.txt

  # not necessary anymore
  #echo "-- convert to table"
  #python2.7 $pythonscripts/uc2otutab.py map.merge.uc > asv_table.merge.txt

  echo "-- create output folder $project.$nowformat; add first files"

  # copy final files into folde
  mkdir -p ../$project.$nowformat
  cp config.txt ../$project.$nowformat/config.txt
  cp zotus.merge.fa ../$project.$nowformat/asvs.merge.fa
  # cp asv_table.merge.txt ../$project.$nowformat/asv_table.merge.txt
  cp asv.tab-csv ../$project.$nowformat/asv_table.merge.txt

  # prepare final project information file
  echo "name;id;own;year;marker;description;participants;doi;repository;accession;ignore" > ../$project.$nowformat/project.csv
  echo "$project;$project;1;$nowformat;$marker;;;;;;" >> ../$project.$nowformat/project.csv

  # prepare final sample information file
  echo "project;name;host;collectionDate;location;country;bioregion;latitude;longitude;tissue;treatment;sampletype;notes" > ../$project.$nowformat/samples.csv
  head -n 1 asv.tab-csv | sed -e "s/#OTU ID[[:space:]]//g" | tr "\t" "\n" | sed "s/^/$project;/"  >> ../$project.$nowformat/samples.csv

fi #end classificationOnly

##### create taxonomy

  echo " "
  echo "===================================="
  echo "Taxonomic classification"
  echo "===================================="

# $vsearch
# old script version: keep to troubleshoot
# refDBs=($(grep "refdb" config.txt | cut -f2 -d"=" | sed 's/\"//g'))
# hieDBs=($(grep "hiedb" config.txt | cut -f2 -d"=" | sed 's/\"//g'))

threshold=$tax_threshold

echo "Classification threshold: $threshold"

# switch here for marker
rm taxonomy.vsearch

if [ "$marker" = "ITS2" ]
  then
    echo ",kingdom,phylum,order,family,genus,species" > taxonomy.vsearch
    echo ",kingdom,phylum,order,family,genus,species" > taxonomy.blast
fi

if [ "$marker" = "COI-5P" ]
  then
    echo ",kingdom,phylum,class,order,family,genus,species" > taxonomy.vsearch
    echo ",kingdom,phylum,class,order,family,subfamily,genus,species, subspecies" > taxonomy.blast
fi

if [ "$marker" = "16S" ]
  then
    echo ",kingdom,phylum,order,family,genus" > taxonomy.vsearch
    echo ",kingdom,phylum,order,family,genus" > taxonomy.blast
fi

countdb=0
cp  zotus.merge.fa zotus.direct.$countdb.uc.nohit.fasta
prevDB=$countdb

for db in "${refDBs[@]}"
  do :
    countdb=$((countdb+1))
    echo "\n\n#### Direct vsearch Classification level: $countdb";
    echo "DB: $db";
    $vsearch --usearch_global zotus.direct.$prevDB.uc.nohit.fasta \
      --db $db \
      --id 0.$threshold \
      --uc zotus.direct.$countdb.uc \
      --fastapairs zotus.direct.$countdb.fasta \
      --strand both \
      --threads $threads 2>  logs/_direct.$countdb.log

    grep "^N[[:space:]]" zotus.direct.$countdb.uc | cut -f 9 > zotus.direct.$countdb.uc.nohit
    $seqfilter zotus.merge.fa --ids zotus.direct.$countdb.uc.nohit --out zotus.direct.$countdb.uc.nohit.fasta
    cut -f 9,10 zotus.direct.$countdb.uc  | grep -v "*" | sed "s/[A-Za-z0-9_-]*;tax=//" >> taxonomy.vsearch
    prevDB=$countdb
  done

echo "\n\n#### Hierarchical vsearch classification";
echo "DB: $hieDBs";

$vsearch --sintax zotus.direct.$countdb.uc.nohit.fasta \
  --db $hieDBs \
  --tabbedout zotus.uc.merge.nohit.sintax \
  --strand plus \
  --sintax_cutoff $sintax_cutoff \
  --threads $threads 2>  logs/_sintax.log

cut -f1,4 zotus.uc.merge.nohit.sintax | sed -E -e "s/\_[0-9]+//g" -e "s/,s:.*$//"  >> taxonomy.vsearch

# v3 idea [TODO]: phylo + spc estimation on sintax results

#cut -f4 zotus.uc.merge.nohit.sintax | head |
#cut -f1,4 zotus.uc.merge.nohit.sintax | head | sed "s/^\(.*\)[[:space:]]\(.*\)$/\2/"

#mkdir -p unresolved

#cut -f4 zotus.uc.merge.nohit.sintax | head  | sort | uniq | sed "s/^.*,//" | sort > unresolved/tmp.taxa

#grep "k:Viridiplantae,p:Streptophyta,c:Magnoliopsida,o:Rosales,f:Moraceae,g:Sorocea" taxonomy.vsearch

## get those with unresolved phylos



# #BLAST local DBs
# countdb=0
# cp  zotus.merge.fa zotus.blast.$countdb.uc.nohit.fasta
# prevDB=$countdb
# touch zotus.blast.hits
# touch taxonomy.blast
#
# for db in "${refDBs[@]}"
#   do :
#     countdb=$((countdb+1))
#     echo "\n\n#### Direct BLAST Classification level: $countdb";
#     echo "DB: $db";
#     #makeblastdb -in $db -parse_seqids -blastdb_version 5 -dbtype nucl
#     blastn  -outfmt '6 qseqid sseqid length pident qcovs' -max_target_seqs 1  -query  zotus.blast.$prevDB.uc.nohit.fasta -subject $db -perc_identity $threshold -qcov_hsp_perc 90 -num_threads $threads > zotus.blast.$countdb.out
#     cut -f1 -d"	" zotus.blast.$countdb.out | cut -f1 >> zotus.blast.hits
#     $sf zotus.merge.fa --ids-exclude --ids  zotus.blast.hits --out zotus.blast.$countdb.uc.nohit.fasta
#     prevDB=$countdb
#     cut -f 1,2 zotus.blast.$countdb.out >> taxonomy.blast
#
#   done

#BLAST web nt instead of sintax here?

echo "-- polishing and copying output files"

sed -i .bak -e "s/c:.*,o:/o:/g" -e "s/,[A-Za-z0-9_-]*;tax=//" -e "s/	/,/"  -e "s/,.*k:/,k:/" taxonomy.vsearch
# sed -i .bak -e "s/c:.*,o:/o:/g" -e "s/[A-Za-z0-9]*;tax=//" -e "s/	/,/" taxonomy.blast

if [ "$classificationOnly" = "16S" ]
  then
    sed -e "s/,.*d:/,d:/"  -e "s/p:Actinobacteria/p:Actinobacteria_phylum/" -e "s/_\d,/,/" taxonomy.vsearch  >  ../$project.$nowformat/taxonomy.vsearch
  else
    sed -e "s/,.*d:/,d:/" taxonomy.vsearch  >  ../$project.$nowformat/taxonomy.vsearch
fi

#sed -i .bak "s/OTUId//" ../$project.$nowformat/asv_table.merge.txt

# without python scripts
sed -i .bak "s/^\#OTU ID//" ../$project.$nowformat/asv_table.merge.txt
sed -i .bak "s/;size=[0-9]*//" ../$project.$nowformat/taxonomy.vsearch

#cp config.txt  ../$project.$nowformat/

echo "-- create R script"

# CREATE R SCRIPT
cat ../_resources/R_template_header.R > ../$project.$nowformat/R_$project.v0.R
echo "# Created: $now" >>../$project.$nowformat/R_$project.v0.R
echo "# Project: $project" >>../$project.$nowformat/R_$project.v0.R
echo "# Marker: $marker" >>../$project.$nowformat/R_$project.v0.R
echo "# For: $datasupplier" >>../$project.$nowformat/R_$project.v0.R

cat ../_resources/R_template_libraries.R >> ../$project.$nowformat/R_$project.v0.R
echo "\n\n# Setting working directory (check path)" >>../$project.$nowformat/R_$project.v0.R
echo "setwd('$(pwd)/../$project.$nowformat')" >>../$project.$nowformat/R_$project.v0.R

echo "\n\n# Custom functions inclusion" >>../$project.$nowformat/R_$project.v0.R
echo "marker=\"$marker\"" >>../$project.$nowformat/R_$project.v0.R
echo "source('$(pwd)/../_resources/metabarcoding_tools_0-1a.R')" >>../$project.$nowformat/R_$project.v0.R

cat ../_resources/R_template_ITS2.R >> ../$project.$nowformat/R_$project.v0.R

cp ../_resources/metabarcoding_tools_0-1a.R ../$project.$nowformat/
mkdir -p ../$project.$nowformat/plots

$vsearch -v > logs/software.version
cd ..

# cleanup if wanted
if [ $compressionCleanup -eq 1 ]
  then
    echo "-- compressing large files"
    sh _compression_cleanup_1.sh $project
fi
