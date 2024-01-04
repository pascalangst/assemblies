##~/software/ccs m64156_230929_124402.subreads.bam m64156_230929_124402.hifi-reads.bam --hifi-kinetics
##
##~/software/jasmine m64156_230929_124402.hifi-reads.bam m64156_230929_124402.5mc.hifi-reads.bam 
##
##minimap2 -ax map-hifi -Y GB-ELK1-1.busco.blast.self.fa <( samtools fastq -TMM,ML m64156_230929_124402.bc2002--bc2002.hifi_reads.fastq.gz ) > GB-ELK1-1.Y.sam
##
##samtools view -@ 5 -Sb GB-ELK1-1.Y.sam | samtools sort -@ 5 -o - > GB-ELK1-1.Y.bam
##rm GB-ELK1-1.Y.sam
##
##~/programs/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
##  --bam GB-ELK1-1.Y.bam \
##  --output-prefix GB-ELK1-1.Y \
##  --model ~/programs/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
##  --threads 8
##  
##  
##  
  
  
  
  
  
  
  
  
~/programs/pbmm2 align GB-ELK1-1.busco.blast.self.fa m64156_230929_124402.bc2002--bc2002.bam GB-ELK1-1.pbmm2.bam --sort -j 8 -J 5

~/programs/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam GB-ELK1-1.pbmm2.bam \
  --output-prefix GB-ELK1-1.pbmm2 \
  --model ~/programs/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 8
  