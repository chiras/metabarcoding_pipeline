#!/bin/bash

###################################################################
# Bioinformatic data processing for metabarcoding
# by Alexander Keller (LMU München)
# email: keller@bio.lmu.de
# if you use this script, please kindly cite this article:
# https://doi.org/10.1098/rstb.2021.0171
###################################################################


githash=$(git log -1 --pretty=format:%h)
gitversion=$(git log --pretty=format:%h | wc -l)
set -o pipefail

if [ -z "$1" ]; then
  echo 'No directory supplied' >&2
  exit 1
fi

project=$(echo "$1" | sed 's:/*$::' | sed 's/^.*\///')
path=$1

now=$(date)
nowformat=$(date +"%y-%m-%d")

mkdir -p $project.$nowformat

# start logging
(
echo "------------------------------------------------------------------"
echo "script: https://github.com/chiras/metabarcoding_pipeline"
echo "version: $githash, revision $gitversion" 
echo "------------------------------------------------------------------"
echo "project: $project"
echo "project path: $path"
echo "results output path: $project.$nowformat"
echo "script log file path: $project.$nowformat/script_$nowformat.log"
echo "------------------------------------------------------------------"

echo "-- moving in data dir"
cd $project

echo "-- loading config"
source config.txt

echo "-- creating directories"
mkdir -p logs
mkdir -p raw
mkdir -p tmp

# Compression tool for intermediate files
if command -v pigz >/dev/null 2>&1; then
    gzip_cmd="pigz -p $threads -1"
else
    gzip_cmd="gzip -1"
fi

if [ $classificationOnly -ne 1 ]
  then
if [ $skip_preprocessing -ne 1 ]
  then

mkdir -p logs/preprocessing
cd raw
echo "-- decompressing raw files"
find . -name '*.gz' -print0 | xargs -0 -I {} -P $threads gunzip {}

# Define the header of the log file
echo "sample_name,project,total_reads,merged_reads,filtered_reads,truncated_reads,filter_strategy,filter_file,final_reads,*sample_name,sample_title,*organism,*collection_date,*env_broad_scale,*env_local_scale,*env_medium,*geo_loc_name,*host,*lat_lon,elev,source_material_id,description" > ../logs/_consolidated_log.csv

# Loop through all files for merging
echo "-- Starting merging and filtering --"
for f in *_R1_*.fastq; do

    # Derive R2 filename
    r="${f/_R1_/_R2_}"

    # Check that R2 exists
    if [[ ! -f "$r" ]]; then
        echo "ERROR: Missing reverse read for $f"
        echo "Expected: $r"
        continue
    fi

    # Derive sample name robustly for both naming schemes
    s="$f"
    s="${s%.fastq}"
    s="${s%.fq}"
    s="${s/_L001_R1_001/}"
    s="${s/_R1_001/}"
    s="${s/_R1_/}"

    total=$(( $(wc -l < "$f") / 4 ))

    echo "===================================="
    echo "Processing sample $s"
    echo "(F: $f R: $r)"
    echo "===================================="

   
    # Stripping last bases for quality
    $vsearch --fastq_filter $f \
        --fastq_trunclen_keep $cutend_fw \
        --fastq_truncee 3 \
        --fastqout $f.sl.fq \
        --threads $threads 2> ../logs/preprocessing/vsearch.rtf.fw.$s.log

    $vsearch --fastq_filter $r \
        --fastq_trunclen_keep $cutend_rv \
        --fastq_truncee 3 \
        --fastqout $r.sl.fq \
        --threads $threads 2> ../logs/preprocessing/vsearch.rtf.rv.$s.log

    # Merging reads
    $vsearch --fastq_mergepairs $f.sl.fq \
          --reverse $r.sl.fq \
          $mergeoptions \
          --fastq_minovlen $merge_minovlen \
          --fastq_maxdiffs $merge_maxdiffs \
          --fastqout $s.merge.fq \
          --fastq_eeout \
          --relabel R1+2-${s}_ \
          --threads $threads  2> ../logs/preprocessing/vsearch.m.$s.log

    merged_reads=$(grep "Merged" ../logs/preprocessing/vsearch.m.$s.log )
    echo "$s : $merged_reads" | tee -a ../logs/_merging.log
    merged_reads=$(echo $merged_reads | grep -o '[0-9]*' | head -n 1)

    # Quality filtering on merged reads
    $vsearch --fastq_filter $s.merge.fq \
      --fastq_stripleft $stripleft \
      --fastq_stripright $stripright \
      --fastq_maxee $fastq_maxee  \
      --fastq_minlen $fastq_minlen  \
      --fastq_maxlen $fastq_maxlen  \
      --fastq_maxns $fastq_maxns \
      --fastaout $s.mergefiltered.fa  \
      --fasta_width 0 \
      --threads $threads 2> ../logs/preprocessing/vsearch.mf.$s.log

    filtered_reads=$(grep "sequences kept" ../logs/preprocessing/vsearch.mf.$s.log )
    echo "$s : $filtered_reads" | tee -a ../logs/_filter.log
    filtered_reads=$(echo $filtered_reads | grep -o '[0-9]*' | head -n 1)

    # Truncation filtering
    $vsearch --fastq_truncee $fastq_truncee \
          --fastq_filter $f \
          --fastq_minlen $fastq_minlen2 \
          --fastaout $s.trunc.fa \
          --relabel R1-${s}_ \
          --threads $threads 2> ../logs/preprocessing/vsearch.tf.$s.log

    trunc_reads=$(grep "sequences kept" ../logs/preprocessing/vsearch.tf.$s.log )
    echo "$s : $trunc_reads" | tee -a ../logs/_truncfilter.log
    trunc_reads=$(echo $trunc_reads | grep -o '[0-9]*' | head -n 1)

    # Apply conditions for selecting the best file
    selected_file=""
    if [[ $filtered_reads -gt 10000 ]]; then
        selected_file="$s.mergefiltered.fa"
        selected_reads=$filtered_reads
        not_selected_reads=$trunc_reads
        strategy="MERGE & Filter"
    elif [[ $filtered_reads -gt 5000 && $trunc_reads -gt 10000 ]]; then
        selected_file="$s.trunc.fa"
        selected_reads=$trunc_reads
        not_selected_reads=$filtered_reads
        strategy="TRUNC & Filter"
    elif [[ $filtered_reads -gt 5000 && $trunc_reads -lt 10000 ]]; then
        selected_file="$s.mergefiltered.fa"
        selected_reads=$filtered_reads
        not_selected_reads=$trunc_reads
        strategy="MERGE & Filter"
    elif [[ $filtered_reads -lt 5000 && $trunc_reads -gt 5000 ]]; then
        selected_file="$s.trunc.fa"
        selected_reads=$trunc_reads
        not_selected_reads=$filtered_reads
        strategy="TRUNC & Filter"
    else
        if [[ $trunc_reads -lt $(( 3 * filtered_reads / 2 )) ]]; then
            selected_file="$s.mergefiltered.fa"
            selected_reads=$filtered_reads
            not_selected_reads=$trunc_reads
            strategy="MERGE & Filter (Low)"
        else
            selected_file="$s.trunc.fa"
            selected_reads=$trunc_reads
            not_selected_reads=$filtered_reads
            strategy="TRUNC & Filter (Low)"
        fi
    fi

    # Output the selected file and copy it as final selection
    cp "$selected_file" "../$s.selection.fa"

    if [[ "$compress_intermediates" -eq 1 ]]; then
        echo "Compressing and moving intermediate files for $s"

        mkdir -p ../tmp/preprocessing

        for tmpfile in \
            "$f.sl.fq" \
            "$r.sl.fq" \
            "$s.merge.fq" \
            "$s.mergefiltered.fa" \
            "$s.trunc.fa"
        do
            if [[ -f "$tmpfile" && "$tmpfile" != "$selected_file" ]]; then

                # Compress if not already compressed
                if [[ ! -f "$tmpfile.gz" ]]; then
                    $gzip_cmd "$tmpfile"
                    tmpfile="$tmpfile.gz"
                else
                    tmpfile="$tmpfile.gz"
                fi

                # Move compressed intermediate to tmp folder
                if [[ -f "$tmpfile" ]]; then
                    mv "$tmpfile" "../tmp/preprocessing/"
                fi
            fi
        done
    fi

    # Log results to consolidated CSV file
    echo "$s,$project,$total,$merged_reads,$filtered_reads,$trunc_reads,$strategy,$selected_file,$selected_reads,$sample_name,$sample_name,metagenome,,,,,,,,,," >> ../logs/_consolidated_log.csv

    # Display selection summary
    echo "Selected $strategy: $selected_file ($selected_reads reads)"
    echo ""
done
cd ..
cat *selection.fa > all.merge.fasta

fi #end skip preprocessing


#     # cleanup
#     echo "-- cleanup"
#     mv *.fastq raw/
#     mv *.fq tmp/
#     mv *.fa tmp/

#   echo "-- removing primer sequences"
#   if [ $skip_primerremoval -ne 1 ] 
#     then
#       python ../_resources/python/remove_primers_2.py \
#         --input all.merge.fasta \
#         --output all.merge.fasta.noprimer.fasta \
#         --marker "$marker" \
#         --max-5p-prefix-stagger "$max_5p_prefix_stagger" \
#         --max-3p-suffix-stagger "$max_3p_suffix_stagger"
#   else
#     echo "Skipping primer removal because skip_primerremoval == 1"
#     #cp all.merge.fasta all.merge.fasta.noprimer.fasta
#   fi

#   echo " "
#   echo "===================================="
#   echo "ASV generation and mapping"
#   echo "===================================="

#   echo "-- derep"

#   derep_input="all.merge.fasta.noprimer.fasta"
#   derep_output="all.merge.derep.fa"
#   derep_uc="all.merge.derep.uc"

#   # Number of FASTA records above which chunked dereplication is used
#   # 2-line FASTA is assumed because primer removal writes fasta_width 0
#   derep_chunk_threshold="${derep_chunk_threshold:-50000000}"

#   # Number of FASTA records per chunk
#   derep_chunk_records="${derep_chunk_records:-10000000}"

#   derep_records=$(grep -c "^>" "$derep_input")

#   echo "Dereplication input records: $derep_records"

#   if [[ "$derep_records" -gt "$derep_chunk_threshold" ]]; then
#       echo "Large input detected. Using chunked dereplication."

#       mkdir -p tmp/derep_chunks tmp/derep_out logs/derep_chunks

#       rm -f tmp/derep_chunks/chunk_* tmp/derep_out/*.derep.fa all.merge.derep.stage1.fa

#       # FASTA-aware split by record count.
#       # This is safe for both wrapped and unwrapped FASTA.
#       awk -v records_per_file="$derep_chunk_records" \
#           -v outdir="tmp/derep_chunks" '
#           BEGIN {
#               file_index = 0
#               record_count = 0
#           }
#           /^>/ {
#               if (record_count % records_per_file == 0) {
#                   if (out) close(out)
#                   file_index++
#                   out = sprintf("%s/chunk_%05d.fa", outdir, file_index)
#               }
#               record_count++
#           }
#           {
#               print > out
#           }
#       ' "$derep_input"

#       total_chunks=$(find tmp/derep_chunks -type f -name "chunk_*.fa" | wc -l | tr -d ' ')
#       chunk_i=0

#       echo "Created $total_chunks dereplication chunks."

#       if [[ "$total_chunks" -eq 0 ]]; then
#           echo "ERROR: splitting produced no chunks."
#           exit 1
#       fi

#       for chunk in tmp/derep_chunks/chunk_*.fa; do
#           chunk_i=$((chunk_i + 1))
#           base=$(basename "$chunk" .fa)

#           echo "-- Chunk dereplication $chunk_i / $total_chunks: $base"

#           "$vsearch" --derep_fulllength "$chunk" \
#             --minuniquesize 1 \
#             --sizein \
#             --sizeout \
#             --fasta_width 0 \
#             --output "tmp/derep_out/${base}.derep.fa" \
#             2> "logs/derep_chunks/${base}.log"

#           if [[ $? -ne 0 || ! -s "tmp/derep_out/${base}.derep.fa" ]]; then
#               echo "ERROR: chunk dereplication failed for $chunk"
#               echo "Check logs/derep_chunks/${base}.log"
#               exit 1
#           fi
#       done


#       echo "-- Concatenating chunk-level dereplicated files"
#       cat tmp/derep_out/*.derep.fa > all.merge.derep.stage1.fa
#       $gzip_cmd tmp/derep_chunks/chunk_*.fa
#       $gzip_cmd tmp/derep_out/*.derep.fa

#       echo "-- Final global dereplication after chunking"

#       "$vsearch" --derep_fulllength all.merge.derep.stage1.fa \
#         --minuniquesize 2 \
#         --sizein \
#         --sizeout \
#         --fasta_width 0 \
#         --uc "$derep_uc" \
#         --output "$derep_output" \
#         2> logs/_derep.log
#       $gzip_cmd all.merge.derep.stage1.fa

#   else
#       echo "Input below chunk threshold. Using standard dereplication."

#       "$vsearch" --derep_fulllength "$derep_input" \
#         --minuniquesize 2 \
#         --sizein \
#         --sizeout \
#         --fasta_width 0 \
#         --uc "$derep_uc" \
#         --output "$derep_output" \
#         2> logs/_derep.log
#   fi

#   if [[ $? -ne 0 || ! -s "$derep_output" ]]; then
#       echo "ERROR: dereplication failed. Check logs/_derep.log"
#       exit 1
#   fi

#   echo "-- denoise"
#   $vsearch --cluster_unoise  all.merge.derep.fa \
#     --sizein --sizeout \
#     --relabel ASV \
#     --unoise_alpha $unoisealpha \
#     --centroids asvs.merge_chim.fa \
#     --threads $threads 2> logs/_unoise.log

#   grep "Clusters" logs/_unoise.log
#   grep "Singletons" logs/_unoise.log

#   echo "-- sorting"
#   $vsearch --sortbysize asvs.merge_chim.fa \
#     --output asvs.merge_sorted.fa \
#     --threads $threads 2>  logs/_sort.log

#   echo "-- chimera removal"
#   $vsearch --uchime3_denovo asvs.merge_sorted.fa \
#     --abskew 16 \
#     --nonchimeras asvs.merge.fa \
#     --threads $threads 2>  logs/_uchime.log

#   grep "Found" logs/_uchime.log

#   if [[ $postcluster -gt 0 ]]; then
#     mv asvs.merge.fa asvs.merge.tmp.fa
#     echo "-- postclustering of ASVs at $postcluster% identity"
#     $vsearch --cluster_fast asvs.merge.tmp.fa \
#               --uc asvs.postcluster.uc \
#               --centroids asvs.merge.fa  \
#               --sizein --sizeout --sizeorder \
#               --threads $threads \
#               --id 0.$postcluster 2>  logs/_postcluster.log
#     clusterin=$(grep -c ">" asvs.merge.tmp.fa)
#     clusterout=$(grep -c ">" asvs.merge.fa)
#     echo "$clusterin ASVs post-clustered to $clusterout pASVs"
#   fi

#   ### create community table
#   echo "-- add barcodelabel"
#   #cat all.merge.fasta.noprimer.fa |  sed "s/^>R1+2-\(.*\)\_\([0-9]*\);/>R1+2-\1_\2;barcodelabel=\1;/g" |  sed "s/^>R1-\([a-zA-Z0-9-]*\)\_\([0-9]*\)/>R1-\1_\2;barcodelabel=\1;/g" > all.merge.bc.fasta
#   python ../_resources/python/add_barcodelabel.py -i all.merge.fasta.noprimer.fasta -o all.merge.bc.fasta

#  echo "-- map data against ASVs"
#  $vsearch --usearch_global all.merge.bc.fasta --db asvs.merge.fa --strand plus --id 0.97 --uc map.merge.uc --otutabout asv.tab-csv --sizeout --threads $threads 2> logs/_mapping.log

  grep "Matching" logs/_mapping.log

  echo "-- create output folder $project.$nowformat; add first files"

  # copy final files into folde
  mkdir -p ../$project.$nowformat
  cp config.txt ../$project.$nowformat/config.txt
  cp asvs.merge.fa ../$project.$nowformat/asvs.merge.fa
  # cp asv_table.merge.txt ../$project.$nowformat/asv_table.merge.txt
  cp asv.tab-csv ../$project.$nowformat/asv_table.merge.txt

  # prepare final project information file
  echo "name;id;own;year;marker;description;participants;doi;repository;accession;ignore" > ../$project.$nowformat/project.csv
  echo "$project;$project;1;$nowformat;$marker;;;;;;" >> ../$project.$nowformat/project.csv

  # prepare final sample information file
  echo "project;name;host;collectionDate;location;country;bioregion;latitude;longitude;tissue;treatment;sampletype;notes" > ../$project.$nowformat/samples.csv
  head -n 1 asv.tab-csv | sed -e "s/#OTU ID[[:space:]]//g" | tr "\t" "\n" | sed "s/^/$project;/"  >> ../$project.$nowformat/samples.csv
  cp logs/_consolidated_log.csv ../$project.$nowformat/samples_metadata.csv
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
if [ "$marker" = "fITS" ]
  then
    echo ",kingdom,phylum,class,order,family,genus,species" > taxonomy.vsearch
    #echo ",kingdom,phylum,order,family,genus,species" > taxonomy.blast
fi

if [ "$marker" = "fITS+16S" ]
  then
    echo ",kingdom,phylum,class,order,family,genus,species" > taxonomy.vsearch
    #echo ",kingdom,phylum,order,family,genus,species" > taxonomy.blast
fi

if [ "$marker" = "ITS2" ]
  then
    echo ",kingdom,phylum,order,family,genus,species" > taxonomy.vsearch
    #echo ",kingdom,phylum,order,family,genus,species" > taxonomy.blast
fi

if [ "$marker" = "COI-5P" ]
  then
    echo ",kingdom,phylum,class,order,family,genus,species" > taxonomy.vsearch
    #echo ",kingdom,phylum,class,order,family,subfamily,genus,species, subspecies" > taxonomy.blast
fi

if [ "$marker" = "16S" ]
  then
    echo ",kingdom,phylum,order,family,genus" > taxonomy.vsearch
    #echo ",kingdom,phylum,order,family,genus" > taxonomy.blast
fi

countdb=0
cp  asvs.merge.fa asvs.direct.$countdb.uc.nohit.fasta
prevDB=$countdb

for db in "${refDBs[@]}"
  do :
    countdb=$((countdb+1))
    echo "-- Direct vsearch Classification level: $countdb";
    echo "DB: $db";
    $vsearch --usearch_global asvs.direct.$prevDB.uc.nohit.fasta \
      --db $db \
      --id 0.$threshold \
      --uc asvs.direct.$countdb.uc \
      --fastapairs asvs.direct.$countdb.fasta \
      --strand both \
      --threads $threads 2>  logs/_direct.$countdb.log

    grep "^N[[:space:]]" asvs.direct.$countdb.uc | cut -f 9 > asvs.direct.$countdb.uc.nohit
    total_xhits=$(wc -l asvs.direct.$countdb.uc | grep -o '[0-9]*' | head -n 1)
    counted_hits=$(grep -c -v "^N[[:space:]]" asvs.direct.$countdb.uc )
    counted_nohits=$(grep -c  "^N[[:space:]]" asvs.direct.$countdb.uc )
    listed_nohits=$(wc -l asvs.direct.$countdb.uc.nohit | grep -o '[0-9]*' | head -n 1)

    echo "Total:                      $total_xhits, thereof $counted_hits hits and $counted_nohits no hits"
    echo "Listed for next iteration:  $listed_nohits"

    python ../_resources/python/subset_fasta.py -i asvs.merge.fa \
      -l asvs.direct.$countdb.uc.nohit \
      -o asvs.direct.$countdb.uc.nohit.fasta
    
    #$seqfilter asvs.merge.fa --ids asvs.direct.$countdb.uc.nohit --out asvs.direct.$countdb.uc.nohit.fasta
    filtered_nohits=$(grep -c  ">" asvs.direct.$countdb.uc.nohit.fasta)
    echo "Filtered:                   $filtered_nohits"

    cut -f 9,10 asvs.direct.$countdb.uc  | grep -v "*" | sed "s/[A-Za-z0-9_-]*;tax=//" >> taxonomy.vsearch
    prevDB=$countdb
  done

echo "-- Hierarchical vsearch classification"

# hieDBs may be supplied either as a single scalar path:
#   hieDBs="/path/to/db.fa"
# or as an indexed array, analogous to refDBs:
#   declare -a hieDBs
#   hieDBs[1]="/path/to/db.level1.fa"
#   hieDBs[2]="/path/to/db.level2.fa"
hieDB_list=()
if declare -p hieDBs 2>/dev/null | grep -q "declare \\-a"; then
    for db in "${hieDBs[@]}"; do
        if [[ -n "$db" ]]; then
            hieDB_list+=("$db")
        fi
    done
else
    if [[ -n "${hieDBs:-}" ]]; then
        hieDB_list+=("$hieDBs")
    fi
fi

if [[ "${#hieDB_list[@]}" -eq 0 ]]; then
    echo "ERROR: no hierarchical reference database supplied via hieDBs"
    exit 1
fi

hier_input="asvs.direct.$countdb.uc.nohit.fasta"
hier_count=0

for hdb in "${hieDB_list[@]}"; do
    hier_count=$((hier_count+1))
    hier_prefix="asvs.hier.$hier_count"

    if [[ ! -s "$hier_input" ]]; then
        echo "No ASVs left for hierarchical classification. Stopping before level $hier_count."
        break
    fi

    input_records=$(grep -c ">" "$hier_input")
    echo "-- Hierarchical vsearch classification level: $hier_count"
    echo "DB: $hdb"
    echo "Input ASVs: $input_records"

    "$vsearch" --sintax "$hier_input" \
      --db "$hdb" \
      --tabbedout "$hier_prefix.sintax" \
      --strand plus \
      --sintax_cutoff "$sintax_cutoff" \
      --threads "$threads" 2> "logs/_sintax.$hier_count.log"

    python ../_resources/python/sintax_overview.py "$hier_prefix.sintax"

    classified_ids="$hier_prefix.classified.ids"
    nohit_ids="$hier_prefix.nohit.ids"
    nohit_fasta="$hier_prefix.nohit.fasta"

    if [ "${use_blast_sintax_combination:-0}" -eq 1 ]; then
        echo "-- LCA BLAST classification level: $hier_count"

        hdb_prefix="${hdb%.*}"

        if [[ ! -f "${hdb_prefix}.ndb" ]]; then
            echo "BLAST database not found. Creating with makeblastdb..."

            if [[ ! -x "$makeblastdb" ]]; then
                echo "ERROR: makeblastdb not found or not executable: $makeblastdb"
                exit 1
            fi

            "$makeblastdb" -in "$hdb" -dbtype nucl -out "$hdb_prefix"
        else
            echo "BLAST database found. Skipping creation."
        fi

        if [[ ! -x "$blastn" ]]; then
            echo "ERROR: blastn not found or not executable: $blastn"
            exit 1
        fi

        "$blastn" -query "$hier_input" -db "$hdb_prefix" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
            -max_target_seqs 50 \
            -num_threads "$threads" \
            -out "$hier_prefix.blast_output.tsv"

        python ../_resources/python/infer_lca.py \
            "$hier_prefix.blast_output.tsv" \
            "$hier_prefix.blast_output.tsv.lca"

        greedy=""
        if [ "${use_blast_sintax_combination_greedy:-0}" -eq 1 ]; then
            greedy="--greedy"
        fi

        python ../_resources/python/combine_taxonomy.py \
            --blast "$hier_prefix.blast_output.tsv.lca" \
            --sintax "$hier_prefix.sintax" \
            --output "$hier_prefix.blast-lca_sintax.out" \
            $greedy

        awk -F '\t' 'NF >= 2 && $2 != "" && $2 != "*" {print $1}' \
            "$hier_prefix.blast-lca_sintax.out" | sort -u > "$classified_ids"

        if [[ -s "$classified_ids" ]]; then
            awk -F '\t' 'NF >= 2 && $2 != "" && $2 != "*" {print $1 "\t" $2}' \
                "$hier_prefix.blast-lca_sintax.out" | \
                sed -E -e "s/_[0-9]+//g" -e "s/,s:.*$//" >> taxonomy.vsearch
        else
            echo "WARNING: no classified IDs detected in $hier_prefix.blast-lca_sintax.out"
            echo "         Falling back to SINTAX-only output for hierarchical level $hier_count."

            cut -f1,4 "$hier_prefix.sintax" | \
                awk -F '\t' 'NF >= 2 && $2 != "" && $2 != "*" {print $0}' | \
                sed -E -e "s/_[0-9]+//g" -e "s/,s:.*$//" >> taxonomy.vsearch

            awk -F '\t' 'NF >= 4 && $4 != "" && $4 != "*" {print $1}' \
                "$hier_prefix.sintax" | sort -u > "$classified_ids"
        fi
    else
        cut -f1,4 "$hier_prefix.sintax" | \
            awk -F '\t' 'NF >= 2 && $2 != "" && $2 != "*" {print $0}' | \
            sed -E -e "s/_[0-9]+//g" -e "s/,s:.*$//" >> taxonomy.vsearch

        awk -F '\t' 'NF >= 4 && $4 != "" && $4 != "*" {print $1}' \
            "$hier_prefix.sintax" | sort -u > "$classified_ids"
    fi

    # Build the input FASTA for the next hierarchical classifier from IDs that
    # remain unclassified at this level.
    cut -f1 "$hier_prefix.sintax" | sort -u > "$hier_prefix.input.ids"
    comm -23 "$hier_prefix.input.ids" "$classified_ids" > "$nohit_ids"

    python ../_resources/python/subset_fasta.py \
      -i "$hier_input" \
      -l "$nohit_ids" \
      -o "$nohit_fasta"

    classified_count=$(wc -l < "$classified_ids" | tr -d ' ')
    nohit_count=$(grep -c ">" "$nohit_fasta" 2>/dev/null || echo 0)

    echo "Classified at hierarchical level $hier_count: $classified_count"
    echo "Remaining for next level:                 $nohit_count"

    hier_input="$nohit_fasta"
done

cp "$hier_input" asvs.hier.final.nohit.fasta


# v3 idea [TODO]: phylo + spc estimation on sintax results

echo "-- polishing and copying output files"

cp taxonomy.vsearch taxonomy.vsearch.bak

python3 ../_resources/python/fix_output_files.py --tax taxonomy.vsearch --asv asv.tab-csv
cp taxonomy.vsearch ../$project.$nowformat/taxonomy.vsearch
cp asv.tab-csv ../$project.$nowformat/asv_table.merge.txt

echo "-- create R script"

# CREATE R SCRIPT
cat ../_resources/R_template_header.R > ../$project.$nowformat/R_$project.v0.R
echo "# Created: $now" >>../$project.$nowformat/R_$project.v0.R
echo "# Project: $project" >>../$project.$nowformat/R_$project.v0.R
echo "# Marker: $marker" >>../$project.$nowformat/R_$project.v0.R
echo "# For: $datasupplier" >>../$project.$nowformat/R_$project.v0.R

cat ../_resources/R_template_libraries.R >> ../$project.$nowformat/R_$project.v0.R
echo "# Setting working directory (check path)" >>../$project.$nowformat/R_$project.v0.R
echo "setwd('$(pwd)/../$project.$nowformat')" >>../$project.$nowformat/R_$project.v0.R

echo "# Custom functions inclusion" >>../$project.$nowformat/R_$project.v0.R
echo "marker=\"$marker\"" >>../$project.$nowformat/R_$project.v0.R
echo "source('./metabarcoding_tools_0-1a.R')" >>../$project.$nowformat/R_$project.v0.R

cat ../_resources/R_template_ITS2.R >> ../$project.$nowformat/R_$project.v0.R

cp ../_resources/metabarcoding_tools_0-1a.R ../$project.$nowformat/
mkdir -p ../$project.$nowformat/plots

if [ $phylogeny -eq 1 ]
  then
    $mafft asvs.merge.fa > asvs.merge.mafft.fa
    sed "s/;size/_size/g" asvs.merge.mafft.fa > asvs.merge.mafft2.fa
    $raxml --msa asvs.merge.mafft2.fa --model GTR+G --msa-format FASTA -threads $threads --prefix asvs
    sed "s/;size/_size/g" ./taxonomy.vsearch > ../$project.$nowformat/taxonomy.vsearch
    sed "s/;size/_size/g" asv.tab-csv > ../$project.$nowformat/asv_table.merge.txt
    cp ./asvs.merge.mafft2.fa.raxml.bestTree ../$project.$nowformat/asvs.tre
fi

$vsearch -v > logs/software.version
cd ..

# cleanup if wanted
if [ $compressionCleanup -eq 1 ]
  then
    echo "-- compressing large files"
    sh _compression_cleanup_1.sh $project
  else 
    echo "NO compression and cleanup applied. Consider calling:"
    echo " bash _compression_cleanup_1.sh $project"
fi

) | tee -a $project.$nowformat/script_$nowformat.log