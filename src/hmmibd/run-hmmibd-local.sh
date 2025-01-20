#!/bin/bash
hmmIBD -i out/hmmIBD-barcodes.tsv -o out/IBD-barcodes -m 150

hmmIBD -i out/hmmIBD-WGS-snps.tsv -o out/IBD-WGS-snps -m 150
