#!/bin/zsh

# for f in $(ls -l | grep "^d" | grep "ITS2" | rev | cut -f1 -d" " | rev); do sh ./_compression_cleanup.sh $f ; done
cd $1

echo " "
echo "------------------------------"
echo "Folder $1"
echo "--- Finalizing with cleanup: "
echo "* Compression raw folder"
tar -zcf raw.tar.gz raw
rm -r raw

echo "* Compression tmp folder"
tar -zcf tmp.tar.gz tmp
rm -r tmp

echo "* Compression logs folder"
tar -zcf logs.tar.gz logs
rm -r logs

echo "* Compression large files"
tar -zcf all.merge.bc.fasta.tar.gz all.merge.bc.fasta
tar -zcf all.merge.derep.uc.tar.gz all.merge.derep.uc
tar -zcf all.merge.fasta.tar.gz all.merge.fasta
tar -zcf all.trunc.fasta.tar.gz all.trunc.fasta
tar -zcf map.merge.uc.tar.gz map.merge.uc

rm all.merge.bc.fasta
rm all.merge.derep.uc
rm all.merge.fasta
rm all.trunc.fasta
rm map.merge.uc
