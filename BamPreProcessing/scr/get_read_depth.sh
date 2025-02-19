#start time
date

#input config
name=H1975_gDNA
echo $name

REF=example_data/Homo_sapiens_assembly38.fasta
BAM=output_bam/$name.sorted.RGadded.marked.realigned.fixed.recal.bam

#output congif
OUTPUTDIR=output_read_depth/
if [ ! -d $OUTPUTDIR ]; then
  mkdir $OUTPUTDIR
fi

#module config
gatkDIR=programs/[GenomeAnalysisTK.jar]
java1_8DIR=programs/[java]

$java1_8DIR -jar $gatkDIR \
 -T DepthOfCoverage \
 -R $REF \
 -I $BAM \
 -o $OUTPUTDIR
# RNA data의 경우 --filter_reads_with_N_cigar option 넣기

#end time
date