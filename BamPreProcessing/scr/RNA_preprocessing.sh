#start time
date

#input config
REF=DNA/example_data/Homo_sapiens_assembly38.fasta
FASTQDIR=DNA/example_data/
FASTQ1=$FASTQDIR/DRR016869_1.tr.fq.gz
FASTQ2=$FASTQDIR/DRR016869_2.tr.fq.gz
dbSNP=/data/references/hg38_gcp/Homo_sapiens_assembly38.dbsnp138.vcf

#output congif
OUTPUTDIR=DNA/output_data/
if [ ! -d $OUTPUTDIR ]; then
  mkdir $OUTPUTDIR
fi
name=H1975_gDNA
PROJECT=test_project
PLATFORM=Illumina
echo $name

#module config
picardDIR=[picard.jar]
gatkDIR=[GenomeAnalysisTK.jar]
java1_8DIR=[java]
#module load picard/2.8.0
#module load java/jdk-1.8u112
#module load bwa/0.7.17

#make log folder
if [ ! -d log ]; then
  mkdir log
fi
if [ ! -d log/tmp ]; then
  mkdir log/tmp
fi


# 1. mapping
bwa mem -t 4 -M $REF $FASTQ1 $FASTQ2 | samtools view -b -o $OUTPUTDIR/$name.bam
echo "bam file done"
date

# 2. sorting
java -XX:+UseSerialGC -Xmx16g -Djava.io.tmpdir=./log/tmp -jar $picardDIR SortSam \
 SO=coordinate \
 INPUT=$OUTPUTDIR/$name.bam \
 OUTPUT=$OUTPUTDIR/$name.sorted.bam \
 VALIDATION_STRINGENCY=LENIENT \
 CREATE_INDEX=true
echo "sorting done"
date

# 3.  Add or replace read groups
java -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $picardDIR AddOrReplaceReadGroups \
 INPUT=$OUTPUTDIR/$name.sorted.bam \
 OUTPUT=$OUTPUTDIR/$name.sorted.RGadded.bam \
 SORT_ORDER=coordinate \
 RGLB=$PROJECT \
 RGPL=$PLATFORM \
 RGPU=$PLATFORM \
 RGSM=$name \
 CREATE_INDEX=true \
 VALIDATION_STRINGENCY=LENIENT
echo "RGadded done"
date

# 3. split N cigar reads
gatk SplitNCigarReads \
 -R $REF \
 -I $OUTPUTDIR/$name.sorted.RGadded.bam \
 -O $OUTPUTDIR/$name.sorted.RGadded.splited.bam
echo "SplitNCigarReads done"
date

# 4. Quality score recalibration. CountCovariates is not available in the new GATK 2.0 version. Instead, use 'BaseRecalibrator' followed by 'PrintReads'
# 4.1) BaseRecalibrator
$java1_8DIR -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $gatkDIR \
 -T BaseRecalibrator \
 -R $REF \
 -I $OUTPUTDIR/$name.sorted.RGadded.splited.bam \
 --knownSites $dbSNP \
 -o $OUTPUTDIR/$name.sorted.RGadded.splited.recal_data.grp \
 --filter_reads_with_N_cigar

# 4.2) PrintReads
$java1_8DIR -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $gatkDIR \
 -T PrintReads \
 -R $REF \
 -I $OUTPUTDIR/$name.sorted.RGadded.splited.bam \
 -BQSR $OUTPUTDIR/$name.sorted.RGadded.splited.recal_data.grp \
 -o $OUTPUTDIR/$name.sorted.RGadded.splited.recal.bam \
 --filter_reads_with_N_cigar
echo "BQSR done"

#end time
date