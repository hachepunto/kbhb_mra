#!/bin/bash
# ARACNE-AP network inference — Basal BRCA: TCGA + METABRIC
# Requires: Java, ARACNe-AP jar (https://github.com/califano-lab/ARACNe-AP)
#
# Inputs  (generados por basal_pre_networks.R):
#   data/tcga_basal_matrix.txt   — expresión Basal TCGA  (genes × muestras, tab-sep)
#   data/mtbrc_basal_matrix.txt  — expresión Basal METABRIC
#   data/TFs_Lambert2018.txt     — lista de TFs Lambert et al. 2018
#
# Outputs (usados por mra_kbhb.R):
#   data/tcga_basal_network.txt  — red consolidada TCGA
#   data/mtbrc_basal_network.txt — red consolidada METABRIC
#
# Runtime: ~2 h por red con 75 cores; ajustar N_CORES y MEM según disponibilidad

set -euo pipefail

ARACNE=/STORAGE/genut/software/ARACNe-AP/dist/aracne.jar
TFFILE=data/TFs_Lambert2018.txt
N_BOOTSTRAPS=100
N_CORES=75
MEM_BOOT=4G      # por proceso bootstrap (75 × 4G = 300 GB en paralelo)
MEM_FULL=32G     # para calculateThreshold y consolidate

# Flags JVM que limitan cada proceso a 1 thread interno.
# Sin esto, cada Java lanza GC + ForkJoinPool propios y el total
# real de cores excede N_CORES × 2-4.
JVM_SINGLE="-XX:+UseSerialGC -XX:ActiveProcessorCount=1 -Djava.util.concurrent.ForkJoinPool.common.parallelism=1"

# ============================================================
# RED 1: TCGA Basal
# ============================================================

EXPFILE=data/tcga_basal_matrix.txt
OUTDIR=data/tcga_aracne_ap_output

mkdir -p "$OUTDIR"

echo "[$(date)] TCGA — calculando umbral MI..."
java -Xmx${MEM_FULL} -jar "$ARACNE" \
    -e "$EXPFILE" -o "$OUTDIR" --tfs "$TFFILE" \
    --pvalue 1E-8 --seed 1 --calculateThreshold

echo "[$(date)] TCGA — ${N_BOOTSTRAPS} bootstraps con ${N_CORES} cores..."
seq 1 "$N_BOOTSTRAPS" | xargs -P "$N_CORES" -I{} \
    java -Xmx${MEM_BOOT} ${JVM_SINGLE} -jar "$ARACNE" \
        -e "$EXPFILE" -o "$OUTDIR" --tfs "$TFFILE" \
        --pvalue 1E-8 --seed {}

echo "[$(date)] TCGA — consolidando..."
java -Xmx${MEM_FULL} -jar "$ARACNE" \
    -e "$EXPFILE" -o "$OUTDIR" --consolidate

cp "$OUTDIR/network.txt" data/tcga_basal_network.txt
echo "[$(date)] TCGA — red guardada en data/tcga_basal_network.txt"

# ============================================================
# RED 2: METABRIC Basal
# ============================================================

EXPFILE=data/mtbrc_basal_matrix.txt
OUTDIR=data/mtbrc_aracne_ap_output

mkdir -p "$OUTDIR"

echo "[$(date)] METABRIC — calculando umbral MI..."
java -Xmx${MEM_FULL} -jar "$ARACNE" \
    -e "$EXPFILE" -o "$OUTDIR" --tfs "$TFFILE" \
    --pvalue 1E-8 --seed 1 --calculateThreshold

echo "[$(date)] METABRIC — ${N_BOOTSTRAPS} bootstraps con ${N_CORES} cores..."
seq 1 "$N_BOOTSTRAPS" | xargs -P "$N_CORES" -I{} \
    java -Xmx${MEM_BOOT} ${JVM_SINGLE} -jar "$ARACNE" \
        -e "$EXPFILE" -o "$OUTDIR" --tfs "$TFFILE" \
        --pvalue 1E-8 --seed {}

echo "[$(date)] METABRIC — consolidando..."
java -Xmx${MEM_FULL} -jar "$ARACNE" \
    -e "$EXPFILE" -o "$OUTDIR" --consolidate

cp "$OUTDIR/network.txt" data/mtbrc_basal_network.txt
echo "[$(date)] METABRIC — red guardada en data/mtbrc_basal_network.txt"

echo "[$(date)] === ARACNE-AP COMPLETO ==="
