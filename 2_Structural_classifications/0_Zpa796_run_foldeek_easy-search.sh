#!/bin/bash


PREFIX="Zpa796_AF2bestmodels"
INPUT_DIR="/Users/dalsasso/Desktop/Posdoc/CAU/structural_analysis/af2/${PREFIX}/all"
OUTPUT_DIR="./easy-search"

mkdir -p "$OUTPUT_DIR"

### PDB annotation
~/foldseek/bin/foldseek easy-search "${INPUT_DIR}"/*.pdb \
    ~/foldseek/databases/pdb "${OUTPUT_DIR}/${PREFIX}.vs.PDB_foldseek.tsv" tmp \
    -s 7.5 --threads 7 --alignment-type 1 --tmscore-threshold 0.5 --format-mode 4 \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob,taxname,theader,tseq


### CATH50 annotation
~/foldseek/bin/foldseek easy-search "${INPUT_DIR}"/*.pdb \
    ~/foldseek/databases/cath50 "${OUTPUT_DIR}/${PREFIX}.vs.CATH50_foldseek.tsv" tmp \
    -s 7.5 --threads 7 --alignment-type 1 --tmscore-threshold 0.5 --format-mode 4 \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob,taxname,theader,tseq


### AFDB Proteome annotation
~/foldseek/bin/foldseek easy-search "${INPUT_DIR}"/*.pdb \
    ~/foldseek/databases/afdb-prot "${OUTPUT_DIR}/${PREFIX}.vs.AFDB-Proteome_foldseek.tsv" tmp \
    -s 7.5 --threads 7 --alignment-type 1 --tmscore-threshold 0.5 --format-mode 4 \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob,taxname,theader,tseq


### AFDB Swiss-Prot annotation
~/foldseek/bin/foldseek easy-search "${INPUT_DIR}"/*.pdb \
    ~/foldseek/databases/afdb-swiss "${OUTPUT_DIR}/${PREFIX}.vs.AFDB-SwissProt_foldseek.tsv" tmp \
    -s 7.5 --threads 7 --alignment-type 1 --tmscore-threshold 0.5 --format-mode 4 \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob,taxname,theader,tseq


### ECOD40 annotation
~/foldseek/bin/foldseek easy-search "${INPUT_DIR}"/*.pdb \
    ~/foldseek/databases/ecod40 "${OUTPUT_DIR}/${PREFIX}.vs.ECOD40_foldseek.tsv" tmp \
    -s 7.5 --threads 7 --alignment-type 1 --tmscore-threshold 0.5 --format-mode 4 \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob,theader,tseq


### SCOPe40 annotation
~/foldseek/bin/foldseek easy-search "${INPUT_DIR}"/*.pdb \
    ~/foldseek/databases/scope40 "${OUTPUT_DIR}/${PREFIX}.vs.SCOPe40_foldseek.tsv" tmp \
    -s 7.5 --threads 7 --alignment-type 1 --tmscore-threshold 0.5 --format-mode 4 \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob,theader,tseq

# Cleanup
rm -r ./tmp
