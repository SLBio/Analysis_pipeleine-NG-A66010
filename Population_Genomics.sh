
#Re-sequencing of 1671 cottons
#(1) Mapped to the NDM13 reference genome using BWA (v.0.7.17-r1188)

bwa mem -k 32 -t 8 ref.fasta sam_1.fq.gz sam_2.fq.gz 

#(2) Generated a gVCF file for each individual according to the best practices using the Genome Analysis Toolkit (GATK, version v4.5.0.0)

gatk --java-options "-Xmx30G" HaplotypeCaller -R ref.fasta -ERC GVCF -I sam.rmdup.bam -O sam.g.vcf.gz

gatk --java-options "-Xmx50G" CombineGVCFs -R ref.fasta \
   -V vcf.list \
   -O sam.gvcf.gz


gatk --java-options "-Xmx50G" GenotypeGVCFs -R ref.fasta \
   -V sam.gvcf.gz \
   -O sam.vcf.gz \
   -G StandardAnnotation


gatk --java-options "-Xmx10G" SelectVariants -select-type SNP -V sam.vcf.gz -O sam.SNP.vcf.gz
gatk --java-options "-Xmx50G" VariantFiltration \
    -V sam.SNP.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O sam.SNP.filter.vcf.gz

#(3) Constructed a neighbor joining (NJ) tree in TreeBest v1.92

treebest nj -b 100 input_file >nj_tree.nwk




