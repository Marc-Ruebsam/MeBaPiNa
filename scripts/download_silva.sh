#!/bin/bash

SILVA_DIR="$1"

SILVA_VERSION="138"
REMOTE_DIR="https://ftp.arb-silva.de/release_${SILVA_VERSION}/Exports"

FASTA_FILENAME="SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva_trunc.fasta"
FASTALIGN_FILENAME="SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta"
SLVMAP_FILENAME="taxmap_slv_ssu_ref_nr_${SILVA_VERSION}.txt"
TAXLIST_FILENAME="tax_slv_ssu_${SILVA_VERSION}.txt"

mkdir -p "$SILVA_DIR"
cd "$SILVA_DIR"

wget "$REMOTE_DIR/${FASTA_FILENAME}.gz"
wget "$REMOTE_DIR/${FASTA_FILENAME}.gz.md5"
if [ "$(md5sum ./${FASTA_FILENAME}.gz)" = "$(cat ${FASTA_FILENAME}.gz.md5)" ]
then rm ${FASTA_FILENAME}.gz.md5
else echo "MD5 filed for ${FASTA_FILENAME}.gz"; exit 1; fi
gunzip "${FASTA_FILENAME}.gz"
echo "rename ${FASTA_FILENAME} reference.fasta"
mv "${FASTA_FILENAME}" "reference_raw.fasta"

wget "$REMOTE_DIR/${FASTALIGN_FILENAME}.gz"
wget "$REMOTE_DIR/${FASTALIGN_FILENAME}.gz.md5"
if [ "$(md5sum ./${FASTALIGN_FILENAME}.gz)" = "$(cat ${FASTALIGN_FILENAME}.gz.md5)" ]
then rm ${FASTALIGN_FILENAME}.gz.md5
else echo "MD5 filed for ${FASTALIGN_FILENAME}.gz"; exit 1; fi
gunzip "${FASTALIGN_FILENAME}.gz"
echo "rename ${FASTALIGN_FILENAME} reference_aligned.fasta"
mv "${FASTALIGN_FILENAME}" "reference_aligned_raw.fasta"

wget "$REMOTE_DIR/taxonomy/${SLVMAP_FILENAME}.gz"
wget "$REMOTE_DIR/taxonomy/${SLVMAP_FILENAME}.gz.md5"
if [ "$(md5sum ${SLVMAP_FILENAME}.gz)" = "$(cat ${SLVMAP_FILENAME}.gz.md5)" ]
then rm ${SLVMAP_FILENAME}.gz.md5
else echo "MD5 filed for ${SLVMAP_FILENAME}.gz"; exit 1; fi
gunzip "${SLVMAP_FILENAME}.gz"
echo "rename ${SLVMAP_FILENAME} slvmap.txt"
mv "${SLVMAP_FILENAME}" "slvmap.txt"

wget "$REMOTE_DIR/taxonomy/${TAXLIST_FILENAME}.gz"
wget "$REMOTE_DIR/taxonomy/${TAXLIST_FILENAME}.gz.md5"
if [ "$(md5sum ${TAXLIST_FILENAME}.gz)" = "$(cat ${TAXLIST_FILENAME}.gz.md5)" ]
then rm ${TAXLIST_FILENAME}.gz.md5
else echo "MD5 filed for ${TAXLIST_FILENAME}.gz"; exit 1; fi
gunzip "${TAXLIST_FILENAME}.gz"
echo "rename ${TAXLIST_FILENAME} taxlist.txt"
mv "${TAXLIST_FILENAME}" "taxlist.txt"
