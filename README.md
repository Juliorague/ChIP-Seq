# ChIP-Seq Analysis of the Transcription Factor HOMEOBOX GENE 1 in *Arabidopsis thaliana*

## Overview
This repository contains the analysis of ChIP-Seq data to identify genomic regions bound by the transcription factor ATH1 in *Arabidopsis thaliana*. The study aims to understand the gene regulation mechanisms mediated by ATH1.

## Authors
- Julio Ramírez Guerrero
- Julián Román Camacho
- Silvestre Ruano Rodríguez
- Manuel Racero de la Rosa
- Rafael Rubio Ramos

## Background
Understanding the genetic mechanisms and environmental signals that regulate plant growth is vital. This study focuses on identifying the binding sites of ATH1 to elucidate its role in gene regulation.

## Objective
The primary objective is to identify the genomic regions where the ATH1 transcription factor binds, providing insights into the gene regulation associated with ATH1 in *Arabidopsis thaliana*.

## Data
The analysis uses ChIP-Seq data available from GEO under accession number GSE157332. The data includes:
- **ChIP Samples**: 
  - ATH1-GFP_chipseq_rep1
  - ATH1-GFP_chipseq_rep2
  - ATH1-GFP_chipseq_rep3
- **Mock Samples**:
  - WT_control_Ler_chipseq_rep1
  - WT_control_Ler_chipseq_rep2
  - WT_control_Ler_chipseq_rep3

## Methodology
### Workflow
1. **Data Preparation**:
    - Quality control of reads using FastQC.
    - Mapping reads to the *Arabidopsis thaliana* reference genome using Bowtie2.
2. **Peak Calling**:
    - Identification of binding sites using MACS2.
3. **Annotation**:
    - Annotation of peaks with genomic features using ChIPseeker.
4. **Motif Analysis**:
    - Identification of DNA motifs within binding sites using HOMER.
5. **Gene Ontology and Pathway Enrichment**:
    - Functional enrichment analysis using clusterProfiler and pathview.

### Scripts
- **Genome and Annotation Download**:
    ```bash
    #!/bin/bash
    # Download reference genome
    wget -O genome.fa.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    gunzip genome.fa.gz
    # Download annotation
    wget -O annotation.gtf.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz
    gunzip annotation.gtf.gz
    # Build genome index
    bowtie2-build genome.fa index
    ```
- **Sample Processing**:
    ```bash
    #!/bin/bash
    SAMPLE_DIR=$1
    SRR=$2
    NAME=$3
    # Download and process fastq files
    fastq-dump --gzip --split-files $SRR
    if [ -f ${SRR}_2.fastq.gz ]; then
       fastqc ${SRR}_1.fastq.gz
       fastqc ${SRR}_2.fastq.gz
       bowtie2 -x ../../genome/index -1 ${SRR}_1.fastq.gz -2 ${SRR}_2.fastq.gz -S $NAME.sam
    else
       fastqc ${SRR}_1.fastq.gz
       bowtie2 -x ../../genome/index -U ${SRR}_1.fastq.gz -S $NAME.sam
    fi
    samtools sort -o $NAME.bam $NAME.sam
    rm $NAME.sam
    rm *.fastq.gz
    samtools index $NAME.bam
    bamCoverage -bs 5 --normalizeUsing CPM --bam $NAME.bam -o $NAME.bw
    ```
- **Peak Calling**:
    ```bash
    #!/bin/bash
    SAMPLE_DIR=$1
    IP=$2
    CTRL=$3
    NAME=$4
    cd $SAMPLE_DIR
    macs2 callpeak -t $IP -c $CTRL -f BAM --outdir . -n $NAME
    ```

## Results
- **Peak Distribution**: Identified 783 peaks where ATH1 binds.
- **Gene Annotation**: Found 582 genes within promoter regions (<= 2kb from TSS).
- **Motif Enrichment**: Highlighted motifs such as TGATTG associated with Homeobox transcription factors.
- **Functional Enrichment**: Genes involved in cell signaling, communication, and oxygen level regulation.

## Visualization
- **Pie Chart**: Distribution of peaks in gene structures.
- **Bar Graph**: Distribution of distances to TSS.
- **UpSet Plot**: Overlap of peaks across genomic regions.
- **MetaGene Plot**: Binding profiles relative to TSS and TTS.
- **IGV Screenshots**: Visual confirmation of binding sites.

## Conclusion
- ATH1 transcription factor significantly binds to the TGATTG motif.
- ATH1 regulates genes involved in essential processes like cell signaling, oxygen regulation, and plant morphology.
- Validation through comparative analysis of wild-type and ATH1 mutant plants.

## References
- Bencivenga, S. et al. (2021). *ARABIDOPSIS THALIANA HOMEOBOX GENE 1 controls plant architecture by locally restricting environmental responses*. Proc Natl Acad Sci U S A.118(17). PMID: 33888582.
