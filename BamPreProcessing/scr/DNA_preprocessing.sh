#start time
date

#input config
REF=example_data/Homo_sapiens_assembly38.fasta
dbSNP=example_data/Homo_sapiens_assembly38.dbsnp138.vcf
FASTQDIR=example_data/
FASTQ1=$FASTQDIR/DRR016869_1.tr.fq.gz
FASTQ2=$FASTQDIR/DRR016869_2.tr.fq.gz

#output congif
OUTPUTDIR=output_bam/
if [ ! -d $OUTPUTDIR ]; then
  mkdir $OUTPUTDIR
fi
name=H1975_gDNA
PROJECT=test_project
PLATFORM=Illumina
echo $name

#module config
picardDIR=programs/[picard.jar]
gatkDIR=programs/[GenomeAnalysisTK.jar]
java1_8DIR=programs/[java]
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

# 3. Marking PCR duplicates
java -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $picardDIR MarkDuplicates \
 INPUT=$OUTPUTDIR/$name.sorted.RGadded.bam \
 OUTPUT=$OUTPUTDIR/$name.sorted.RGadded.marked.bam \
 METRICS_FILE=metrics \
 CREATE_INDEX=true \
 VALIDATION_STRINGENCY=LENIENT
echo "mark duplicate done"
date

# 4. Local realignment around indels
# [step1] To create a table of possible indels
$java1_8DIR -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $gatkDIR \
 -T RealignerTargetCreator \
 -R $REF \
 -o $OUTPUTDIR/$name.sorted.RGadded.marked.bam.list \
 -I $OUTPUTDIR/$name.sorted.RGadded.marked.bam \
 -dt NONE \
 -nt 8

# [step2] To realign reads around indels targets
$java1_8DIR -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $gatkDIR \
 -I $OUTPUTDIR/$name.sorted.RGadded.marked.bam \
 -R $REF \
 -T IndelRealigner \
 -targetIntervals $OUTPUTDIR/$name.sorted.RGadded.marked.bam.list \
 -o $OUTPUTDIR/$name.sorted.RGadded.marked.realigned.bam \
 --maxReadsForRealignment 300000
echo "local realignment done"
date

# 5. The mate information must be fixed
java -XX:+UseSerialGC -Xmx4g -Djava.io.tmpdir=./log/tmp -jar $picardDIR FixMateInformation \
 INPUT=$OUTPUTDIR/$name.sorted.RGadded.marked.realigned.bam \
 OUTPUT=$OUTPUTDIR/$name.sorted.RGadded.marked.realigned.fixed.bam \
 SO=coordinate \
 VALIDATION_STRINGENCY=LENIENT \
 CREATE_INDEX=true
echo "mate fixing done"
date

# 6. Quality score recalibration. CountCovariates is not available in the new GATK 2.0 version. Instead, use 'BaseRecalibrator' followed by 'PrintReads'
# 6.1) BaseRecalibrator
$java1_8DIR -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $gatkDIR \
 -T BaseRecalibrator \
 -R $REF \
 -I $OUTPUTDIR/$name.sorted.RGadded.marked.realigned.fixed.bam \
 --knownSites $dbSNP \
 -o $OUTPUTDIR/$name.sorted.RGadded.marked.realigned.fixed.bam.recal_data.grp

# 6.2) PrintReads
$java1_8DIR -XX:+UseSerialGC -Xmx36g -Djava.io.tmpdir=./log/tmp -jar $gatkDIR \
 -T PrintReads \
 -R $REF \
 -I $OUTPUTDIR/$name.sorted.RGadded.marked.realigned.fixed.bam \
 -BQSR $OUTPUTDIR/$name.sorted.RGadded.marked.realigned.fixed.bam.recal_data.grp \
 -o $OUTPUTDIR/$name.sorted.RGadded.marked.realigned.fixed.recal.bam
echo "BQSR done"

#end time
date