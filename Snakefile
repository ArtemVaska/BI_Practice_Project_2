URL_ref="http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz"
URL_control_1="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz"
URL_control_2="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz"
URL_control_3="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz"


rule sample_download_from_SRA:  # download roommate's sequencing results
	output:
		"raw_data/roommate.fastq.gz"
	shell:
		"wget -O {output} {URL_ref}"  # -O allows you to save results to a file


rule data_unpack:  # unpacking sequencing results
	input:
		"raw_data/{sample}.fastq.gz"
	output:
		"raw_data/{sample}.fastq"
	shell:
		"gunzip -c {input} > {output}"  # -c prints the result to stdout, and using > we intercept it into the desired file, the original archive remains unchanged


rule reference_download_from_GenBank:  # loading the influenza hemagglutinin gene reference
	output:
		"raw_data/ref/reference.fasta"
	shell:
		"efetch -db nucleotide -id KF848938.1 -format fasta > {output}"


rule control_1_download_from_SRA:  # loading the 1st technical replicate of the control sample
	output:
		"raw_data/control1.fastq.gz"
	shell:
		"wget -O {output} {URL_control_1}"


rule control_2_download_from_SRA:  # loading the 2nd technical replicate of the control sample
	output:
		"raw_data/control2.fastq.gz"
	shell:
		"wget -O {output} {URL_control_2}"


rule control_3_download_from_SRA:  # loading the 3rd technical replicate of the control sample
	output:
		"raw_data/control3.fastq.gz"
	shell:
		"wget -O {output} {URL_control_3}"


rule bwa_index:  # genome indexing -- you get a lot of index files that will be needed for alignment
	input:
		"raw_data/ref/{reference}.fasta"
	log:
		"logs/bwa_index.{reference}.log"
	output:
		"raw_data/ref/{reference}.fasta.amb",
		"raw_data/ref/{reference}.fasta.ann",
		"raw_data/ref/{reference}.fasta.bwt",
		"raw_data/ref/{reference}.fasta.pac",
		"raw_data/ref/{reference}.fasta.sa"
	shell:
		"bwa index {input} &> {log}"


rule bwa_align:  # alignment of reads to reference using bwa and conversion to binary format using samtools
	input:
		"raw_data/ref/{reference}.fasta.amb",
		"raw_data/ref/{reference}.fasta.ann",
		"raw_data/ref/{reference}.fasta.bwt",
		"raw_data/ref/{reference}.fasta.pac",
		"raw_data/ref/{reference}.fasta.sa",
		ref="raw_data/ref/{reference}.fasta",
		reads="raw_data/{sample}.fastq"
	log:
		"logs/bwa_align.{reference}.{sample}.log"
	output:
		temporary("results/bwa/{reference}.{sample}.unsorted.bam")
	threads: 8
	shell:
		"bwa mem -t {threads} {input.ref} {input.reads} 2> {log} | samtools view -S -b > {output}"


rule bam_sort:  # sorting the resulting alignments by chromosome and coordinate within it
	input:
		rules.bwa_align.output
	log:
		"logs/bam_sort.{reference}.{sample}.log"
	output:
		protected("results/bwa/{reference}.{sample}.sorted.bam")
	threads: 8
	shell:
		"samtools sort --threads {threads} {input} > {output}"


rule samtools_index:  # indexing for the ability to open in a GUI utility
	input:
		rules.bam_sort.output
	output:
		"results/bwa/{reference}.{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"


rule mpileup_create:  # creates an mpileup file, which tabulates the number of bases that match or don’t match the reference
	input:
		rules.samtools_index.output,  # mpileup requires a sorted, indexed bam file
		ref="raw_data/ref/{reference}.fasta",
		align=rules.bam_sort.output
	output:
		"results/variants/{reference}.{sample}.mpileup"
	shell:
		"samtools mpileup -d 0 -f {input.ref} {input.align} > {output}"  # since our variants may be quite rare, we set depth limit to maximum with -d 0 flag


rule varscan_freq:  # look for common variants with VarScan
	input:
		rules.mpileup_create.output
	output:
		"results/variants/{reference}.{sample}_freq.vcf"
	shell:
		"java -jar tools/VarScan.v2.4.6.jar mpileup2snp {input} --min-var-freq 0.95 --variants --output-vcf 1 > {output}"


rule varscan_rare:  # look for rare variants with VarScan
	input:
		rules.mpileup_create.output
	output:
		"results/variants/{reference}.{sample}_rare.vcf"
	shell:
		"java -jar tools/VarScan.v2.4.6.jar mpileup2snp {input} --min-var-freq 0.001 --variants --output-vcf 1 > {output}"


rule awk_vcf:  # parsing vcf files using awk tool
	input:
		"results/variants/{name}.vcf"
	output:
		"results/variants/{name}_parsed.txt"
	shell:
		"""cat {input} | awk '{{FS="\t"; OFS="\t"}} NR>24 {{split($10, line, /:/); print $2, $4, $5, substr(line[7], 1, length(line[7])-1)}}' > {output}"""


