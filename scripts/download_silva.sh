#!/bin/bash

KRAKEN2_DB_NAME="$1"

SILVA_VERSION="138"
REMOTE_DIR="https://ftp.arb-silva.de/release_${SILVA_VERSION}/Exports"

FASTA_FILENAME="SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta"
NCBIMAP_FILENAME="taxmap_embl-ebi_ena_ssu_ref_nr99_${SILVA_VERSION}.txt"
SLVMAP_FILENAME="taxmap_slv_ssu_ref_nr_${SILVA_VERSION}.txt"
TAXLIST_FILENAME="tax_slv_ssu_${SILVA_VERSION}.txt"

mkdir -p "$KRAKEN2_DB_NAME"
cd "$KRAKEN2_DB_NAME"

wget "$REMOTE_DIR/${FASTA_FILENAME}.gz"
gunzip "${FASTA_FILENAME}.gz"
mv "${FASTA_FILENAME}" "reference.fasta"
wget "$REMOTE_DIR/taxonomy/${NCBIMAP_FILENAME}.gz"
gunzip "${NCBIMAP_FILENAME}.gz"
mv "${NCBIMAP_FILENAME}" "ncbimap.txt"
wget "$REMOTE_DIR/taxonomy/${SLVMAP_FILENAME}.gz"
gunzip "${SLVMAP_FILENAME}.gz"
mv "${SLVMAP_FILENAME}" "slvmap.txt"
wget "$REMOTE_DIR/taxonomy/${TAXLIST_FILENAME}.gz"
gunzip "${TAXLIST_FILENAME}.gz"
mv "${TAXLIST_FILENAME}" "taxlist.txt"

mkdir -p kraken2/species_tmp/library \
kraken2/species_tmp/taxonomy \
kraken2/genus_tmp/library \
kraken2/genus_tmp/taxonomy \
krona/species \
krona/genus
sed -e '/^>/!y/U/T/' "reference.fasta" > "kraken2/species_tmp/library/library.fna"
cp "kraken2/species_tmp/library/library.fna" "kraken2/genus_tmp/library/library.fna"
