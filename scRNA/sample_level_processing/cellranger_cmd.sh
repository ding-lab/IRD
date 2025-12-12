cellranger=/diskmnt/Software/cellranger-6.0.1/cellranger
reference=/diskmnt/Software/cellranger-6.0.1/References/refdata-gex-GRCh38-2020-A

sample_id=BM2
FASTQ=/diskmnt/primary/MMY10/xfer.genome.wustl.edu/gxfer1/82059335740947/
SAMPLE=H_YX-BM2-Bc1_1-BM2-Bc1_1

chemistry=SC3Pv3

nohup ${cellranger} count \
        --id=${sample_id}\
        --fastqs=${FASTQ}\
	--sample=${SAMPLE}\
        --transcriptome=${reference}\
        --chemistry=${chemistry} \
        --jobmode=local \
        --localcores=12 \
        --localmem=67 &
