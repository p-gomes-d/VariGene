#!/bin/bash

echo ""
echo ""
echo "[[Bem-vind@ ao VariGene]]"
sleep 1.0
echo ""
echo "[[Amostra a analisar: AMOSTRA_A]]"
sleep 1.0
echo "[[Recolher o respectivo ficheiro fastq..]]"
samtools fastq ./files/AMOSTRA_A.bam > AMOSTRA_A.fq 2>&1
echo "[[..sucesso]]"
sleep 1.0
echo "[[Filtrar o ficheiro fastq..]]"
fastp -i AMOSTRA_A.fq -o AMOSTRA_A_cleaned.fq > fastp_stats.txt 2>&1
echo "[[..sucesso]]"
sleep 1.0
echo "[[Descarregar o genoma humano de referência GRCh37 (GenBank Assembly GCA_000001405.1)..]]"
wget hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gzip -d hg19.fa.gz
echo "[[..sucesso]]"
sleep 1.0
echo "[[Indexar o genoma de referência (aguarde durante alguns minutos)..]]"
bwa index -p indexed_hg19 hg19.fa
echo "[[..sucesso]]"
sleep 1.0
echo "[[Alinhar a AMOSTRA_A com o genoma de referência e produzir o ficheiro bam..]]"
bwa mem -M indexed_hg19 AMOSTRA_A_cleaned.fq > AMOSTRA_A.sam
samtools view -Sb AMOSTRA_A.sam > AMOSTRA_A.bam
samtools sort AMOSTRA_A.bam -o AMOSTRA_A_sorted.bam
echo "[[..sucesso]]"
sleep 1.0
echo "[[Indexar o ficheiro bam..]]"
samtools index AMOSTRA_A_sorted.bam
echo "[[..sucesso]]"
sleep 1.0
samtools flagstat AMOSTRA_A_sorted.bam > bwa_stats.txt
samtools faidx hg19.fa
echo "[[Identificar as variantes genéticas..]]"
samtools mpileup -uf hg19.fa AMOSTRA_A_sorted.bam | bcftools call --ploidy GRCh37 -mv -Ov -o AMOSTRA_A.vcf
echo "[[..sucesso]]"
sleep 1.0
echo "[[Filtrar variantes de qualidade (mínimo cobertura de read de 20 e mínimo de qualidade de 20 (probabilidade de erro 0.01 ou inferior))..]]"
bcftools filter -i 'DP>=20 && QUAL>=20' AMOSTRA_A.vcf > AMOSTRA_A_cleaned.vcf
echo "[[..sucesso]]"
sleep 1.0
echo "[[Descarregar a versão mais recente da base de dados ClinVar..]]"
wget ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
gzip -d clinvar.vcf.gz
echo "[[..sucesso]]"
sleep 1.0
echo "[[Preparar a base de dados..]]"
cat <(awk 'NR>=1 && NR<=42 {print}' clinvar.vcf) <(awk 'NR>'$(grep -n "#CHROM" clinvar.vcf | cut -d ":" -f 1)' {print "chr" $0}' clinvar.vcf) > clinvar_withchr.vcf
bgzip -c clinvar_withchr.vcf > clinvar_withchr.vcf.gz
echo "[[..sucesso]]"
sleep 1.0
echo "[[Indexar a base de dados..]]"
bcftools index clinvar_withchr.vcf.gz
echo "[[..sucesso]]"
sleep 1.0
echo "[[Indexar o ficheiro vcf..]]"
bgzip -c AMOSTRA_A_cleaned.vcf > AMOSTRA_A_cleaned.vcf.gz
bcftools index AMOSTRA_A_cleaned.vcf.gz
echo "[[..sucesso]]"
sleep 1.0
echo "[[Anotar o ficheiro vcf através da base de dados ClinVar..]]"
bcftools annotate -a clinvar_withchr.vcf.gz -c CHROM,POS,REF,ALT,INFO/ALLELEID,INFO/CLNDN,INFO/CLNHGVS,INFO/CLNSIG,INFO/CLNVC,INFO/GENEINFO,INFO/MC AMOSTRA_A_cleaned.vcf.gz -Ov -o AMOSTRA_A_annotated.vcf
echo "[[..sucesso]]"
sleep 1.0
echo "[[Gerar o relatório da análise genómica]]"
bash report.sh
echo ""
echo ""
sleep 1.0
echo "[[Obrigad@ por utilizar o VariGene]]"
echo ""
echo ""
