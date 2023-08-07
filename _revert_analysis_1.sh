#!/bin/zsh
if [ -z "$1" ]; then
  echo 'No directory supplied' >&2
  exit 1
fi

cd $1

rm all*
rm asv*
rm -r logs
rm map*
rm tax*
rm -r tmp
rm zotus*
mv raw/* ./
