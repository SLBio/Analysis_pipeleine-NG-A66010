##1. Repetitive sequences identification

#(1) Homology-based prediction
# RepeatMasker (v4.0.5)
RepeatMasker -a -nolow -no_is -norna -parallel 40 -lib ./data/RepeatMasker.lib -s genome.fasta >genome.fa.1.log 2>genome.fa.2.log

#RepeatProteinMask (v4.0.5)
RepeatProteinMask -noLowSimple -pvalue 0.0001 -engine ncbi genome.fasta

#(2) Ab initio prediction

#Using RepeatModeler (v1.0.8) to generate a de novo repeat library (repeatmodeler.library)
BuildDatabase -engine ncbi -name genome genome.fasta
RepeatModeler -database genome -engine ncbi -threads 15 >&run.out

#Using LTR FINDER (v1.0.7) to generate a de novo repeat library (ltr.library)
ltr_finder -C -w 2 -s ./data/Sc-tRNAs.fa genome.fasta >./ltr.library 2>./ltr.log

#Using RepeatScout (v1.0.5) to generate a de novo repeat library (RepeatScout.library)
build_lmer_table -sequence genome.fasta -freq ./genome.freq
RepeatScout -sequence ./genome.fasta -freq ./genome.freq -output ./genome.out
perl ./scripts/filter-stage-1.prl ./genome.out >./RepeatScout.library

#Using PILER (v3.3.0) to generate a de novo repeat library (RepeatScout.library)
sh ./script/Piler.sh genome.fa 20 ./Piler.out

#Using RepeatProteinMask (v4.0.5) generate Repetitive sequences using 4 de novo repeat libraries

cat repeatmodeler.library ltr.library RepeatScout.library RepeatScout.library >denovo.lib
perl ./script/changeAa.pl denovo.lib >denovo.lib.convert
perl ./script/uclust.pl denovo.lib.convert >denovo.lib.convert.unredundance.fa
RepeatMasker -a -nolow -no_is -norna -parallel 40 -lib denovo.lib.convert.unredundance.fa -s genome.fasta >genome.fa.1.log 2>genome.fa.2.log

#(3) Using Identifying tandem repeats sequences using TRF
trf genome.fasta 2 7 7 80 10 50 2000 -d -h



##2. Gene structure

#(1) Homology-based prediction. 

blastall -p tblastn -e 1e-05 -F T -m 8 -d genome.fasta -i query.pep.fasta -o query.pep.blast
perl ./script/solar.pl -a prot2genome2 -z -f m8 query.pep.blast >query.pep.blast.solar
perl ./script/struct_solar_filter.pl query.pep.blast.solar 0.25
perl ./script/struct_solar_to_table.pl query.pep.blast.solar.filter 
perl ./script/struct_genomic_cluster.pl -overlap_percent 0.5 query.pep.blast.solar.filter.table > query.pep.blast.solar.filter.table.nonredundance
perl ./script/fishInWinter.pl -bf table -ff table  query.pep.blast.solar.filter.table.nonredundance query.pep.blast.solar.filter >query.pep.blast.solar.filter.nr
./wise2.4.1/genewise -trev -genesf -gff -sum each.query.gene.fa each.genome.region.fa >each.genome.region.genewise


#(2) Predict genes using RNA-seq data

#predicted the isoform structure (gene model)
./bowtie2-2.2.5/bin/bowtie2-build -q -f genome.fasta genome
./tophat-2.0.13/bin/tophat -p 6 --max-intron-length 500000 -m 2 --library-type fr-unstranded  -o ./outdir genome tissue_1.fq.gz tissue_2.fq.gz
./Cufflinks-2.1.1/cufflinks  -I 500000 -p 1 --library-type fr-unstranded -L CUFF -o ./outdir accepted_hits.bam
perl ./script/Cufflink_gtf2gff3.pl transcripts.gtf >transcripts.gff3

# Assembled transcripts from RNA-seq data using Trinity (v2.1.1)

./trinity-2.1.1/bin/Trinity --seqType fq --left tissue_1.fq.gz --right tissue_2.fq.gz --CPU 20 --max_memory 200G --normalize_reads --full_cleanup --min_glue 2 --min_kmer_cov 2 --KMER_SIZE 25 --output trinity.out


#Assembled transcripts were aligned to the assembled genome using the software PASA to generate the training set
perl ./script/pipeline_pasa.pl --genome_seq genome.fasta  --round 1 --cpu 10 --sql_name genome --trinity_seq trinity.out.fasta --outdir pasa1

perl ./script/training_gff3Togff.pl pasa1.gff3 >pasa1.gff
perl ./script/training_grep_complete.pl pasa1.pep pasa1.complete.pep
perl ./script/training_blast_database.pl --swissprot --evalue 1e-05 --cutf 20 --cpu 20 pasa1.complete.pep
./blast-2.2.26/bin/blastall -b 5  -F F  -p blastp -e 1e-05  -d uniprot_sprot.fasta -a 10 -i pasa1.complete.pep -o pasa1.swissprot.blast
perl script/blast_parser.pl ./pasa1.swissprot.blast > ./pasa1.swissprot.blast.tab;
mv pasa1.swissprot.blast.tab pasa1.complete.pep.nr.blast.tab
awk '$13>=95' pasa1.complete.pep.nr.blast.tab >pasa1.complete.pep.nr.blast.tab.score95
perl ./script/training_blast_stat.pl pasa1.complete.pep.nr.blast.tab.score95 >pasa1.complete.pep.nr.blast.tab.ratio
awk '$4>=0.8' pasa1.complete.pep.nr.blast.tab.ratio >pasa1.complete.pep.nr.blast.tab.ratio.cvg0.8
perl ./script/fastaDeal.pl -attr id:len pasa1.complete.pep >pasa1.complete.pep.len
perl ./script/training_find_nosame_with_head.pl pasa1.complete.pep.nr.blast.tab.ratio pasa1.complete.pep.len >pasa1.complete.pep.nr.blast.tab.ratio.noalign
awk '$2>=1000' pasa1.complete.pep.nr.blast.tab.ratio.noalign >pasa1.complete.pep.nr.blast.tab.ratio.noalign.longORF
cat pasa1.complete.pep.nr.blast.tab.ratio.cvg0.8 pasa1.complete.pep.nr.blast.tab.ratio.noalign.longORF >pasa1.train
perl ./script/training_get.gff.pl pasa1.train pasa1.gff >pasa1.train.gff
perl ./script/training_perfect_gene4pasa.pl --start 10 --stop 10 genome.fasta pasa1.train.gff
perl ./script/training_clustergff.pl pasa1.train.gff
perl ./script/training_trainset_uniq.pl pasa1.train.gff.nr.gff
perl ./script/training_trainset_random.pl pasa1.train.gff.nr.gff.uniq.gff 1000


#(3) Ab initio-based prediction predict using the soft-masked genome using Augustus (v2.5.5), GlimmerHMM (v3.0.1) and SNAP (Semi-HMM-based Nucleic Acid Parser)

#  Augustus (v2.5.5)
perl ./augustus/scripts/autoAugTrain.pl --genome=genome.masked.fasta --trainingset=PASA4training/pasa1.train.gff.nr.gff.uniq.gff.random.gff --optrounds=0 --species=pasa1
./augustus/bin/augustus --species=pasa1 --AUGUSTUS_CONFIG_PATH=training/Augustus/conf --extrinsicCfgFile=training/Augustus/conf/extrinsic/extrinsic.cfg --uniqueGeneId=true --noInFrameStop=true --gff3=on --genemodel=complete --strand=both  genome.masked.fasta >genome.augustus


# GlimmerHMM (v3.0.1)
perl ./script/training_auto_train_glimmerhmmm.pl pasa1.train.gff.nr.gff.uniq.gff.random.gff genome.masked.fasta
./script/glimmerhmm genome.masked.fasta -d training/GlimmerHMM/pasa1 -f -g >genome.gff

# SNAP (Semi-HMM-based Nucleic Acid Parser)
perl ./script/training_train_snap_new.pl --name genome.masked.fasta pasa1 pasa1.train.gff.nr.gff.uniq.gff.random.gff
./script/snap-2013-11-29/snap -gff training/SNAP/pasa1.hmm genome.masked.fasta >genome.snap

#GeneID (v1.4)
geneid -P ./data/homo_sapiens.param -v -G -p geneid genome.masked.fasta >genome.Geneid

#GeneScan (v1.0)
./script/genscan ./data/HumanIso.smat genome.masked.fasta >genome.genscan

#(4) EVidenceModeler (v1.1.1)  was used to integrate all predictions

#To configure a weight file (weights.txt) for running NextDenovo
PROTEIN GeneWise        100
TRANSCRIPT      CUFF    100
ABINITIO_PREDICTION     Augustus        10
ABINITIO_PREDICTION     genscan 1
ABINITIO_PREDICTION     GeneID  1
ABINITIO_PREDICTION     GlimmerHMM      1
ABINITIO_PREDICTION     SNAP    1
OTHER_PREDICTION        transdecoder    100

# Running EVidenceModeler (v1.1.1) 
We used the parameters recommended by the software to run the EVidenceModeler.
Please see https://github.com/EVidenceModeler/EVidenceModeler/wiki#running-evm
We generated the file "evm.gff3".

#(5) PASA2 was used to predict untranslated regions and alternative splicing variations

perl ./script/pipeline_pasa.pl --round 2 --trinity_seq pasa1/genome.assemblies.fasta.transdecoder.cds --genome_seq genome.fasta  --sql_name genome --mysqllib pasa1/mysql_bin/data/  --evm evm.gff3 --outdir 08.pasa 

