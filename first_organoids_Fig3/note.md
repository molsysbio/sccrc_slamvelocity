# How to run this:
1. download the data
2. change file path names accordingly
3. run demux or demux_interactive (you need to demux yourself)
4. run pp_annot
5. run Fig3C_dyn_SLAM_analysis
6. the rest can be run in any order


# cellSNP and vireo
To use cellSNP and vireo you can use this code snippet we used:
```
export BASEDIR=/fast/scratch/users/peidlis_c/sodar_patient_organoid_data

export WDIR=${BASEDIR}/NB_AS_C
cd ${WDIR}

BARCODES=${WDIR}/20200328_cellranger_count/filtered_feature_bc_matrix/barcodes.tsv.gz
if test -f "$BARCODES"; then
    echo "$BARCODES exist"
else
    echo "$BARCODES does not exist, unpacking..."
    tar -xzvf ${WDIR}/20200328_cellranger_count/filtered_feature_bc_matrix.tgz -C ${WDIR}/20200328_cellranger_count/
fi

mkdir -p ${WDIR}/snps/
cellSNP -s ${WDIR}/20200328_cellranger_count/possorted_genome_bam.bam \
	-b ${WDIR}/20200328_cellranger_count/filtered_feature_bc_matrix/barcodes.tsv.gz \
	-Osnps -R /fast/work/users/peidlis_c/utils/single_cell_rna_seq/sample_barcoding/SNPs/genome/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
	-O ${WDIR}/snps/ -p 20 --minMAF 0.1 --minCOUNT 20

mkdir -p ${WDIR}/demuxed/
vireo -c ${WDIR}/snps/ -N 3 -o ${WDIR}/demuxed/
```

# HTO demux with seurat
To demultiplex using Seurat's HTOdemux we refer to the official vignette:
https://satijalab.org/seurat/v3.0/hashing_vignette.html
We followed this as well and used default values for this.
The hashtag <-> identitys are provided in the code.
