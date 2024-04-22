rule trim_se:
	input:
		lambda wildcards: fastq_dir + "/" + samples.Fq1[samples.Name == wildcards.sampleName]
	output:
		trim_dir + "/{sampleName}_se.trim.fq.gz"
	message:
		"Trimming... [{wildcards.sampleName}]"
	log:
		trim_dir + "/{sampleName}.trim.log"
	shell:
		"""
		module load Cutlery

		cutadapt {adapterSeq} {opt_cutadapt} --trim-n \
		-o ${{TMPDIR}}/__temp__.$$.1.fq.gz {input} 2>&1 | tee {log}
		mv ${{TMPDIR}}/__temp__.$$.1.fq.gz {output}
		"""

def get_fastq(sampleName):
	## single-end
	if doTrim:
		return trim_dir + "/{sampleName}_se.trim.fq.gz"
	else:
		return fastq_dir + "/" + samples.Fq1[samples.Name == sampleName].tolist()[0]

rule align_se_rRNA:
	input:
		fq1=lambda wildcards: get_fastq(wildcards.sampleName)
	output:
		summary = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.log",
		bam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.bam.bai",
		unaligned = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/unaligned_pre.fq.gz",
	message:
		"Aligning... [{wildcards.sampleName}]"
	threads:
		cluster["align_se"]["cpu"]
	params:
		desDir = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align"
	shell:
		"""
		module purge
		module load hisat2/2.2.1
		module load samtools/1.18.0

		mkdir -p {params.desDir}

		hisat-3n --un-gz {output.unaligned} {soft_clip} {RNA_splice} --no-unal --summary-file {output.summary} \
		-x {rRNA_index} \
		-q -U {input.fq1} \
		--base-change C,T --repeat --directional-mapping | samtools sort -o {output.bam}

		samtools index {output.bam}
		"""

rule get_align_stat_rRNA:
	input:
		summary = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.log"
	output:
		alignStats = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/alignStats.tsv"
	message:
		"Getting alignment stats... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3

		rBis.alignStat_perSample.py -s {input.summary} -o {output.alignStats} -n {wildcards.sampleName}
		"""

rule post_process_bam_rRNA:
	input:
		fq1 = lambda wildcards: get_fastq(wildcards.sampleName),
		unalignedFq = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/unaligned_pre.fq.gz",
		bam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.bam.bai"
	output:
		newUnaligned = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/unaligned.fq.gz",
		plusBam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.bam.bai",
		plusUniqBam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.uniq.bam",
		plusUniqBamBai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.uniq.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.minus.bam.bai"
	message:
		"Post-processing bam file from rRNA alignment... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		module load seqtk/1.3

		rBis.postProcessAlign.sh \
		-o {sample_dir}/{wildcards.sampleName}/{rRNA_dir}/Align \
		{input.bam} {input.fq1} {input.unalignedFq}

		"""

rule get_read_stats_rRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.bam.bai",
		plusBam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.minus.bam.bai",
		unalignedFQ = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/unaligned.fq.gz",
	output:
		readStats = sample_dir + "/{sampleName}/" + rRNA_dir + "/readStats.tsv"
	message:
		"Counting reads used in mpileup for rRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load samtools/1.18.0
		rBis.getReadStats.sh -o {output.readStats} {input.bam} {input.plusBam} {input.minusBam} {wildcards.sampleName} {input.unalignedFQ}
		"""


rule run_samtools_mpileup_all_rRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.bam.bai"
	output:
		table = sample_dir + "/{sampleName}/" + rRNA_dir + "/mpileup_out_raw_all.tsv"
	message:
		"Running samtools mpileup for rRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {rRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}
		"""

rule run_samtools_mpileup_uniq_rRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.uniq.bam.bai"
	output:
		table = sample_dir + "/{sampleName}/" + rRNA_dir + "/mpileup_out_raw_uniq.tsv"
	message:
		"Running samtools mpileup uniq for rRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {rRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""

rule call_m5C_candidates_rRNA:
	input:
		tableAll = sample_dir + "/{sampleName}/" + rRNA_dir + "/mpileup_out_raw_all.tsv",
		tableUniq = sample_dir + "/{sampleName}/" + rRNA_dir + "/mpileup_out_raw_uniq.tsv"
	output:
		cov = sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv",
		stats = sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv",
		statsSimplified = sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv"
	params:
		cov = cov_thresh,
		desDir = sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir
	message:
		"Calling m5C candidates for rRNA... [{wildcards.sampleName}]"
	shell:
		"""

		module purge
		module load anaconda3
		source activate BisKit

		mkdir -p {params.desDir}

		rBis.mpileup_getStrand.py -c {input.tableAll} -o $TMPDIR/{wildcards.sampleName}_mpileup_all_rRNA.tsv
		rBis.mpileup_getStrand.py -c {input.tableUniq} -o $TMPDIR/{wildcards.sampleName}_mpileup_uniq_rRNA.tsv

		rBis.mpileup_mergeTables.py \
		-a $TMPDIR/{wildcards.sampleName}_mpileup_all_rRNA.tsv \
		-u $TMPDIR/{wildcards.sampleName}_mpileup_uniq_rRNA.tsv \
		-o $TMPDIR/{wildcards.sampleName}_mpileup_final_rRNA.tsv \

		rBis.mpileup_call_m5C.py \
		-v {params.cov} -c $TMPDIR/{wildcards.sampleName}_mpileup_final_rRNA.tsv \
		-o {sample_dir}/{wildcards.sampleName}/{rRNA_dir}/call_{coverage_filtered_dir} \
		-m rRNA -s {sig_thresh_call} -t {sig_type_call} -l {lookup_table} -n {min_C_count} \
		-e {wildcards.sampleName} -u {species} -a {rRNA_fa} -r {genome_fa}

		rm $TMPDIR/{wildcards.sampleName}_mpileup_all_rRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_uniq_rRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_final_rRNA.tsv

		"""

rule draw_distribution_plots_rRNA:
	input:
		cov = sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	output:
		plots = expand(sample_dir + "/{{sampleName}}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/{pre}.png", pre=["log2_cov", "methRate", "log2_cov_"+sig_type_call, "methRate_"+sig_type_call])
	message:
		"Drawing distribution plots for rRNA... [{wildcards.sampleName}]"
	shell:
		"""

		module load python3/3.6.3

		rBis.draw_dist_plots.py \
		-c {input.cov} -o {sample_dir}/{wildcards.sampleName}/{rRNA_dir}/call_{coverage_filtered_dir} \
		-n {wildcards.sampleName} \
		-s {sig_thresh_call} \
		-t {sig_type_call}
		"""


rule align_se_tRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/unaligned.fq.gz"
	output:
		summary = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.log",
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.bam.bai",
		unaligned = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/unaligned_pre.fq.gz"
	message:
		"Aligning tRNA... [{wildcards.sampleName}]"
	threads:
		cluster["align_se"]["cpu"]
	params:
		desDir = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align"
	shell:
		"""
		module purge
		module load hisat2/2.2.1
		module load samtools/1.18.0

		mkdir -p {params.desDir}

		hisat-3n --un-gz {output.unaligned} -k 500 {soft_clip} {RNA_splice} --no-unal --summary-file {output.summary} \
		-x {tRNA_index} \
		-q -U {input.fq1} \
		--base-change C,T --repeat --directional-mapping | samtools sort -o {output.bam}

		samtools index {output.bam}
		"""

rule get_align_stat_tRNA:
	input:
		summary = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.log"
	output:
		alignStats = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/alignStats.tsv"
	message:
		"Getting alignment stats... [{wildcards.sampleName}]"
	shell:
		"""

		module load python3/3.6.3

		rBis.alignStat_perSample.py -s {input.summary} -o {output.alignStats} -n {wildcards.sampleName}
		"""

rule post_process_bam_tRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/unaligned.fq.gz",
		unalignedFq = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/unaligned_pre.fq.gz",
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.bam.bai",
	output:
		newUnaligned = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/unaligned.fq.gz",
		plusBam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam.bai",
		plusUniqBam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam",
		plusUniqBamBai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.minus.bam.bai",
	message:
		"Post-processing bam file from tRNA alignment... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		module load seqtk/1.3

		rBis.postProcessAlign.sh \
		-o {sample_dir}/{wildcards.sampleName}/{tRNA_dir}/Align \
		{input.bam} {input.fq1} {input.unalignedFq}

		"""

rule feature_count_all_tRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + tRNA_dir + "/featureCount_all.txt"
	message:
		"Feature count all tRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3

		rBis.count_chromosomes.py -b {input.bam} > {output.cnt}
		"""

rule feature_count_uniq_tRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam.bai"
	output:
		cnt = sample_dir + "/{sampleName}/" + tRNA_dir + "/featureCount_uniq.txt"
	message:
		"Feature count uniq tRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3

		rBis.count_chromosomes.py -b {input.bam} > {output.cnt}
		"""

rule get_read_stats_tRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.bam.bai",
		plusBam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.minus.bam.bai",
		unalignedFQ = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/unaligned.fq.gz",
	output:
		readStats = sample_dir + "/{sampleName}/" + tRNA_dir + "/readStats.tsv"
	message:
		"Counting reads used in mpileup for tRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load samtools/1.18.0
		rBis.getReadStats.sh -o {output.readStats} {input.bam} {input.plusBam} {input.minusBam} {wildcards.sampleName} {input.unalignedFQ}
		"""


rule run_samtools_mpileup_all_tRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.bam.bai"
	output:
		table = sample_dir + "/{sampleName}/" + tRNA_dir + "/mpileup_out_raw_all.tsv"
	message:
		"Running samtools mpileup for tRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {tRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}
		"""

rule run_samtools_mpileup_uniq_tRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + tRNA_dir + "/mpileup_out_raw_uniq.tsv"
	message:
		"Running samtools mpileup uniq for tRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {tRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""

rule call_m5C_candidates_tRNA:
	input:
		tableAll = sample_dir + "/{sampleName}/" + tRNA_dir + "/mpileup_out_raw_all.tsv",
		tableUniq = sample_dir + "/{sampleName}/" + tRNA_dir + "/mpileup_out_raw_uniq.tsv"
	output:
		cov = sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv",
		stats = sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv",
		statsSimplified = sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv"
	params:
		cov = cov_thresh,
		desDir = sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir
	message:
		"Calling m5C candidates for tRNA... [{wildcards.sampleName}]"

	shell:
		"""

		module purge
		module load anaconda3
		source activate BisKit

		rBis.mpileup_getStrand.py -c {input.tableAll} -o $TMPDIR/{wildcards.sampleName}_mpileup_all_tRNA.tsv
		rBis.mpileup_getStrand.py -c {input.tableUniq} -o $TMPDIR/{wildcards.sampleName}_mpileup_uniq_tRNA.tsv

		rBis.mpileup_mergeTables.py \
		-a $TMPDIR/{wildcards.sampleName}_mpileup_all_tRNA.tsv \
		-u $TMPDIR/{wildcards.sampleName}_mpileup_uniq_tRNA.tsv \
		-o $TMPDIR/{wildcards.sampleName}_mpileup_final_tRNA.tsv \

		rBis.mpileup_call_m5C.py \
		-v {params.cov} -c $TMPDIR/{wildcards.sampleName}_mpileup_final_tRNA.tsv \
		-o {sample_dir}/{wildcards.sampleName}/{tRNA_dir}/call_{coverage_filtered_dir} \
		-m tRNA -s {sig_thresh_call} -t {sig_type_call} -l {lookup_table} -n {min_C_count} \
		-e {wildcards.sampleName} -u {species} -r {genome_fa}

		rm $TMPDIR/{wildcards.sampleName}_mpileup_all_tRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_uniq_tRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_final_tRNA.tsv

		"""

rule draw_distribution_plots_tRNA:
	input:
		cov = sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	output:
		plots = expand(sample_dir + "/{{sampleName}}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/{pre}.png", pre=["log2_cov", "methRate", "log2_cov_"+sig_type_call, "methRate_"+sig_type_call])
	message:
		"Drawing distribution plots for tRNA... [{wildcards.sampleName}]"

	shell:
		"""

		module load python3/3.6.3

		rBis.draw_dist_plots.py \
		-c {input.cov} -o {sample_dir}/{wildcards.sampleName}/{tRNA_dir}/call_{coverage_filtered_dir} \
		-n {wildcards.sampleName} \
		-s {sig_thresh_call} \
		-t {sig_type_call}
		"""


rule align_se_miRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/unaligned.fq.gz"
	output:
		summary = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.log",
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.bam.bai",
		unaligned = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/unaligned_pre.fq.gz"
	message:
		"Aligning miRNA... [{wildcards.sampleName}]"
	threads:
		cluster["align_se"]["cpu"]
	params:
		desDir = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align"
	shell:
		"""
		module purge
		module load hisat2/2.2.1
		module load samtools/1.18.0

		mkdir -p {params.desDir}

		hisat-3n --un-gz {output.unaligned} {soft_clip} {RNA_splice} --no-unal --summary-file {output.summary} \
		-x {miRNA_hairpin_index} \
		-q -U {input.fq1} \
		--base-change C,T --repeat --directional-mapping | samtools sort -o {output.bam}

		samtools index {output.bam}
		"""

rule get_align_stat_miRNA:
	input:
		summary = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.log"
	output:
		alignStats = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/alignStats.tsv"
	message:
		"Getting alignment stats... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.alignStat_perSample.py -s {input.summary} -o {output.alignStats} -n {wildcards.sampleName}
		"""

rule post_process_bam_miRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/unaligned.fq.gz",
		unalignedFq = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/unaligned_pre.fq.gz",
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.bam.bai"
	output:
		newUnaligned = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/unaligned.fq.gz",
		plusBam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam.bai",
		plusUniqBam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam",
		plusUniqBamBai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.minus.bam.bai"
	message:
		"Post-processing bam file from miRNA alignment... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		module load seqtk/1.3

		rBis.postProcessAlign.sh \
		-o {sample_dir}/{wildcards.sampleName}/{miRNA_dir}/Align \
		{input.bam} {input.fq1} {input.unalignedFq}

		"""

rule feature_count_all_miRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam.bai"
	output:
		cnt = sample_dir + "/{sampleName}/" + miRNA_dir + "/featureCount_all.txt",
		summary = sample_dir + "/{sampleName}/" + miRNA_dir + "/featureCount_all.txt.summary",
		cntlog = sample_dir + "/{sampleName}/" + miRNA_dir + "/featureCount_all.log"
	message:
		"Feature count all miRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load subread/2.0.3
		rBis.featureCountsSample.sh -o {sample_dir}/{wildcards.sampleName}/{miRNA_dir}/featureCount_all -g {miRNA_gtf} -s "-M -O -T 1 -t mature -g gene_id --fracOverlap 0.8" {input.bam}
		"""

rule feature_count_uniq_miRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam.bai"
	output:
		cnt = sample_dir + "/{sampleName}/" + miRNA_dir + "/featureCount_uniq.txt",
		summary = sample_dir + "/{sampleName}/" + miRNA_dir + "/featureCount_uniq.txt.summary",
		cntlog = sample_dir + "/{sampleName}/" + miRNA_dir + "/featureCount_uniq.log"
	message:
		"Feature count uniq miRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load subread/2.0.3
		rBis.featureCountsSample.sh -o {sample_dir}/{wildcards.sampleName}/{miRNA_dir}/featureCount_uniq -g {miRNA_gtf} -s "-O -T 1 -t mature -g gene_id --fracOverlap 0.8" {input.bam}
		"""

rule get_read_stats_miRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.bam.bai",
		plusBam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.minus.bam.bai",
		unalignedFQ = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/unaligned.fq.gz",
	output:
		readStats = sample_dir + "/{sampleName}/" + miRNA_dir + "/readStats.tsv"
	message:
		"Counting reads used in mpileup for miRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load samtools/1.18.0
		rBis.getReadStats.sh -o {output.readStats} {input.bam} {input.plusBam} {input.minusBam} {wildcards.sampleName} {input.unalignedFQ}
		"""

rule run_samtools_mpileup_all_miRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.bam.bai"
	output:
		table = sample_dir + "/{sampleName}/" + miRNA_dir + "/mpileup_out_raw_all.tsv"
	message:
		"Running samtools mpileup for miRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {miRNA_hairpin_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}
		"""

rule run_samtools_mpileup_uniq_miRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + miRNA_dir + "/mpileup_out_raw_uniq.tsv"
	message:
		"Running samtools mpileup uniq for miRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {miRNA_hairpin_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""

rule call_m5C_candidates_miRNA:
	input:
		tableAll = sample_dir + "/{sampleName}/" + miRNA_dir + "/mpileup_out_raw_all.tsv",
		tableUniq = sample_dir + "/{sampleName}/" + miRNA_dir + "/mpileup_out_raw_uniq.tsv"
	output:
		cov = sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv",
		stats = sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv",
		statsSimplified = sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv"
	params:
		cov = cov_thresh,
		desDir = sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir
	message:
		"Calling m5C candidates for miRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load anaconda3
		source activate BisKit

		mkdir -p {params.desDir}

		rBis.mpileup_getStrand.py -c {input.tableAll} -o $TMPDIR/{wildcards.sampleName}_mpileup_all_miRNA.tsv
		rBis.mpileup_getStrand.py -c {input.tableUniq} -o $TMPDIR/{wildcards.sampleName}_mpileup_uniq_miRNA.tsv

		rBis.mpileup_mergeTables.py \
		-a $TMPDIR/{wildcards.sampleName}_mpileup_all_miRNA.tsv \
		-u $TMPDIR/{wildcards.sampleName}_mpileup_uniq_miRNA.tsv \
		-o $TMPDIR/{wildcards.sampleName}_mpileup_final_miRNA.tsv \

		rBis.mpileup_call_m5C.py \
		-v {params.cov} -c $TMPDIR/{wildcards.sampleName}_mpileup_final_miRNA.tsv \
		-o {sample_dir}/{wildcards.sampleName}/{miRNA_dir}/call_{coverage_filtered_dir} \
		-m miRNA -s {sig_thresh_call} -t {sig_type_call} -l {lookup_table} -n {min_C_count} \
		-e {wildcards.sampleName} -u {species} -r {genome_fa} -g {miRNA_gtf}

		rm $TMPDIR/{wildcards.sampleName}_mpileup_all_miRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_uniq_miRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_final_miRNA.tsv

		"""

rule draw_distribution_plots_miRNA:
	input:
		cov = sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	output:
		plots = expand(sample_dir + "/{{sampleName}}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/{pre}.png", pre=["log2_cov", "methRate", "log2_cov_"+sig_type_call, "methRate_"+sig_type_call])
	message:
		"Drawing distribution plots for miRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module load python3/3.6.3
		rBis.draw_dist_plots.py \
		-c {input.cov} -o {sample_dir}/{wildcards.sampleName}/{miRNA_dir}/call_{coverage_filtered_dir} \
		-n {wildcards.sampleName} \
		-s {sig_thresh_call} \
		-t {sig_type_call}
		"""

rule align_se_piRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/unaligned.fq.gz"
	output:
		summary = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.log",
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.bam.bai",
		unaligned = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/unaligned_pre.fq.gz"
	message:
		"Aligning piRNA... [{wildcards.sampleName}]"
	threads:
		cluster["align_se"]["cpu"]
	params:
		desDir = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align"
	shell:
		"""
		module purge
		module load hisat2/2.2.1
		module load samtools/1.18.0
		
		mkdir -p {params.desDir}

		hisat-3n --un-gz {output.unaligned} {soft_clip} {RNA_splice} --no-unal --summary-file {output.summary} \
		-x {piRNA_index} \
		-q -U {input.fq1} \
		--base-change C,T --repeat --directional-mapping | samtools sort -o {output.bam}

		samtools index {output.bam}
		"""

rule get_align_stat_piRNA:
	input:
		summary = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.log"
	output:
		alignStats = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/alignStats.tsv"
	message:
		"Getting alignment stats... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.alignStat_perSample.py -s {input.summary} -o {output.alignStats} -n {wildcards.sampleName}
		"""

rule post_process_bam_piRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/unaligned.fq.gz",
		unalignedFq = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/unaligned_pre.fq.gz",
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.bam.bai",
	output:
		newUnaligned = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/unaligned.fq.gz",
		plusBam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam.bai",
		plusUniqBam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam",
		plusUniqBamBai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.minus.bam.bai",
	message:
		"Post-processing bam file from piRNA alignment... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		module load seqtk/1.3

		rBis.postProcessAlign.sh \
		-o {sample_dir}/{wildcards.sampleName}/{piRNA_dir}/Align \
		{input.bam} {input.fq1} {input.unalignedFq}

		"""

rule feature_count_all_piRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + piRNA_dir + "/featureCount_all.txt"
	message:
		"Feature count all piRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.count_chromosomes.py -b {input.bam} > {output.cnt}
		"""

rule feature_count_uniq_piRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + piRNA_dir + "/featureCount_uniq.txt"
	message:
		"Feature count uniq piRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.count_chromosomes.py -b {input.bam} > {output.cnt}
		"""

rule get_read_stats_piRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.bam.bai",
		plusBam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.minus.bam.bai",
		unalignedFQ = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/unaligned.fq.gz",

	output:
		readStats = sample_dir + "/{sampleName}/" + piRNA_dir + "/readStats.tsv"
	message:
		"Counting reads used in mpileup for piRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load samtools/1.18.0
		rBis.getReadStats.sh -o {output.readStats} {input.bam} {input.plusBam} {input.minusBam} {wildcards.sampleName} {input.unalignedFQ}
		"""

rule run_samtools_mpileup_all_piRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + piRNA_dir + "/mpileup_out_raw_all.tsv"
	message:
		"Running samtools mpileup for piRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		samtools mpileup -B -f {piRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}
		"""

rule run_samtools_mpileup_uniq_piRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + piRNA_dir + "/mpileup_out_raw_uniq.tsv"
	message:
		"Running samtools mpileup uniq for piRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		samtools mpileup -B -f {piRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""

rule call_m5C_candidates_piRNA:
	input:
		tableAll = sample_dir + "/{sampleName}/" + piRNA_dir + "/mpileup_out_raw_all.tsv",
		tableUniq = sample_dir + "/{sampleName}/" + piRNA_dir + "/mpileup_out_raw_uniq.tsv"
	output:
		cov = sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv",
		stats = sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv",
		statsSimplified = sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv"
	params:
		cov = cov_thresh,
		desDir = sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir
	message:
		"Calling m5C candidates for piRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load anaconda3
		source activate BisKit

		mkdir -p {params.desDir}

		rBis.mpileup_getStrand.py -c {input.tableAll} -o $TMPDIR/{wildcards.sampleName}_mpileup_all_piRNA.tsv
		rBis.mpileup_getStrand.py -c {input.tableUniq} -o $TMPDIR/{wildcards.sampleName}_mpileup_uniq_piRNA.tsv

		rBis.mpileup_mergeTables.py \
		-a $TMPDIR/{wildcards.sampleName}_mpileup_all_piRNA.tsv \
		-u $TMPDIR/{wildcards.sampleName}_mpileup_uniq_piRNA.tsv \
		-o $TMPDIR/{wildcards.sampleName}_mpileup_final_piRNA.tsv \

		rBis.mpileup_call_m5C.py \
		-v {params.cov} -c $TMPDIR/{wildcards.sampleName}_mpileup_final_piRNA.tsv \
		-o {sample_dir}/{wildcards.sampleName}/{piRNA_dir}/call_{coverage_filtered_dir} \
		-m piRNA -s {sig_thresh_call} -t {sig_type_call} -l {lookup_table} -n {min_C_count} \
		-e {wildcards.sampleName}  -u {species} -r {genome_fa}

		rm $TMPDIR/{wildcards.sampleName}_mpileup_all_piRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_uniq_piRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_final_piRNA.tsv

		"""

rule draw_distribution_plots_piRNA:
	input:
		cov = sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	output:
		plots = expand(sample_dir + "/{{sampleName}}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/{pre}.png", pre=["log2_cov", "methRate", "log2_cov_"+sig_type_call, "methRate_"+sig_type_call])
	message:
		"Drawing distribution plots for piRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module load python3/3.6.3
		rBis.draw_dist_plots.py \
		-c {input.cov} -o {sample_dir}/{wildcards.sampleName}/{piRNA_dir}/call_{coverage_filtered_dir} \
		-n {wildcards.sampleName} \
		-s {sig_thresh_call} \
		-t {sig_type_call}
		"""


rule align_se_genome:
	input:
		fq1 = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/unaligned.fq.gz"
	output:
		summary = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.log",
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai",
		unaligned = sample_dir + "/{sampleName}/" + genome_dir + "/Align/unaligned.fq.gz"
	message:
		"Aligning genome... [{wildcards.sampleName}]"
	threads:
		cluster["align_se"]["cpu"]
	params:
		desDir = sample_dir + "/{sampleName}/" + genome_dir + "/Align"
	shell:
		"""
		module purge
		module load hisat2/2.2.1
		module load samtools/1.18.0

		mkdir -p {params.desDir}

		hisat-3n --un-gz {output.unaligned} {soft_clip} {genome_splice} --no-unal --summary-file {output.summary} \
		-x {genome_index} \
		-q -U {input.fq1} \
		--base-change C,T --repeat --directional-mapping | samtools sort -o {output.bam}

		samtools index {output.bam}

		"""

rule get_align_stat_genome:
	input:
		summary = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.log"
	output:
		alignStats = sample_dir + "/{sampleName}/" + genome_dir + "/Align/alignStats.tsv"
	message:
		"Getting alignment stats... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.alignStat_perSample.py -s {input.summary} -o {output.alignStats} -n {wildcards.sampleName}
		"""

rule post_process_bam_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai",
	output:
		uniqBam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.uniq.bam",
		uniqBamBai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.uniq.bam.bai",
		plusBam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.minus.bam.bai"
	message:
		"Post-processing bam file from genome alignment... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		## split to forward strand
		echo "making align.plus.bam."
		samtools view -b -F 16 {input.bam} > {output.plusBam}
		samtools index {output.plusBam}

		## split to minus strand
		echo "making align.minus.bam."
		samtools view -b -f 16 {input.bam} > {output.minusBam}
		samtools index {output.minusBam}

		## get align.uniq.bam
		echo "making align.uniq.bam."
		samtools view -H {input.bam} > {sample_dir}/{wildcards.sampleName}/{genome_dir}/Align/header.txt
		samtools view {input.bam} | grep -w "NH:i:1" > $TMPDIR/{wildcards.sampleName}.genome.tmp
		cat {sample_dir}/{wildcards.sampleName}/{genome_dir}/Align/header.txt $TMPDIR/{wildcards.sampleName}.genome.tmp > $TMPDIR/{wildcards.sampleName}.genome.tmp.sam
		samtools view -bS $TMPDIR/{wildcards.sampleName}.genome.tmp.sam > {output.uniqBam}
		samtools index {output.uniqBam}

		rm {sample_dir}/{wildcards.sampleName}/{genome_dir}/Align/header.txt
		rm $TMPDIR/{wildcards.sampleName}.genome.tmp
		rm $TMPDIR/{wildcards.sampleName}.genome.tmp.sam

		"""

rule feature_count_all_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + genome_dir + "/featureCount_all.txt",
		summary = sample_dir + "/{sampleName}/" + genome_dir + "/featureCount_all.txt.summary",
		cntlog = sample_dir + "/{sampleName}/" + genome_dir + "/featureCount_all.log"
	message:
		"Feature count all genome... [{wildcards.sampleName}]"
	shell:
		"""
		module load subread/2.0.3
		rBis.featureCountsSample.sh -o {sample_dir}/{wildcards.sampleName}/{genome_dir}/featureCount_all -g {genome_bed} -s "-M -O -T 1 -t exon -g gene_name --fracOverlap 0.8" {input.bam}
		"""

rule feature_count_uniq_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.uniq.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + genome_dir + "/featureCount_uniq.txt",
		summary = sample_dir + "/{sampleName}/" + genome_dir + "/featureCount_uniq.txt.summary",
		cntlog = sample_dir + "/{sampleName}/" + genome_dir + "/featureCount_uniq.log"
	message:
		"Feature count uniq genome... [{wildcards.sampleName}]"
	shell:
		"""
		module load subread/2.0.3
		rBis.featureCountsSample.sh -o {sample_dir}/{wildcards.sampleName}/{genome_dir}/featureCount_uniq -g {genome_bed} -s "-O -T 1 -t exon -g gene_name --fracOverlap 0.8" {input.bam}
		"""

rule get_read_stats_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai",
		plusBam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.minus.bam.bai",
		unalignedFQ = sample_dir + "/{sampleName}/" + genome_dir + "/Align/unaligned.fq.gz",

	output:
		readStats = sample_dir + "/{sampleName}/" + genome_dir + "/readStats.tsv"
	message:
		"Counting reads used in mpileup for genome... [{wildcards.sampleName}]"
	shell:
		"""
		module load samtools/1.18.0
		rBis.getReadStats.sh -o {output.readStats} {input.bam} {input.plusBam} {input.minusBam} {wildcards.sampleName} {input.unalignedFQ}
		"""

rule run_samtools_mpileup_all_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai"
	output:
		table = sample_dir + "/{sampleName}/" + genome_dir + "/mpileup_out_raw_all.tsv"
	message:
		"Running samtools mpileup for genome... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		samtools mpileup -B -f {genome_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""

rule run_samtools_mpileup_uniq_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.uniq.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + genome_dir + "/mpileup_out_raw_uniq.tsv"
	message:
		"Running samtools mpileup uniq for genome... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		samtools mpileup -B -f {genome_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""


rule call_m5C_candidates_genome:
	input:
		tableAll = sample_dir + "/{sampleName}/" + genome_dir + "/mpileup_out_raw_all.tsv",
		tableUniq = sample_dir + "/{sampleName}/" + genome_dir + "/mpileup_out_raw_uniq.tsv"
	output:
		cov = sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir + "/call.tsv",
		stats = sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv",
		statsSimplified = sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv"
	params:
		cov = cov_thresh,
		desDir = sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir
	message:
		"Calling m5C candidates for genome... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load anaconda3
		source activate BisKit

		mkdir -p {params.desDir}

		rBis.mpileup_getStrand.py -c {input.tableAll} -o $TMPDIR/{wildcards.sampleName}_mpileup_all_genome.tsv
		rBis.mpileup_getStrand.py -c {input.tableUniq} -o $TMPDIR/{wildcards.sampleName}_mpileup_uniq_genome.tsv

		rBis.mpileup_mergeTables.py \
		-a $TMPDIR/{wildcards.sampleName}_mpileup_all_genome.tsv \
		-u $TMPDIR/{wildcards.sampleName}_mpileup_uniq_genome.tsv \
		-o $TMPDIR/{wildcards.sampleName}_mpileup_final_genome.tsv \

		gawk -F'\t' '{{ if (NR!=1) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5 }}' $TMPDIR/{wildcards.sampleName}_mpileup_final_genome.tsv > $TMPDIR/{wildcards.sampleName}_{genome_dir}_1.seqContext_input.bed
		bedtools intersect -loj -a $TMPDIR/{wildcards.sampleName}_{genome_dir}_1.seqContext_input.bed -b {genome_bed} > $TMPDIR/{wildcards.sampleName}_{genome_dir}_2.annotated.bed

		rBis.mpileup_call_m5C.py \
		-v {params.cov} -c $TMPDIR/{wildcards.sampleName}_{genome_dir}_2.annotated.bed \
		-y $TMPDIR/{wildcards.sampleName}_mpileup_final_genome.tsv \
		-o {sample_dir}/{wildcards.sampleName}/{genome_dir}/call_{coverage_filtered_dir} \
		-m genome -s {sig_thresh_call} -t {sig_type_call} -l {lookup_table} -n {min_C_count} \
		-e {wildcards.sampleName} -u {species} -r {genome_fa}

		rm $TMPDIR/{wildcards.sampleName}_{genome_dir}_1.seqContext_input.bed
		rm $TMPDIR/{wildcards.sampleName}_{genome_dir}_2.annotated.bed
		rm $TMPDIR/{wildcards.sampleName}_mpileup_all_genome.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_uniq_genome.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_final_genome.tsv

		"""

rule draw_distribution_plots_genome:
	input:
		cov = sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	output:
		plots = expand(sample_dir + "/{{sampleName}}/" + genome_dir + "/call_" + coverage_filtered_dir + "/{pre}.png", pre=["log2_cov", "methRate", "log2_cov_"+sig_type_call, "methRate_"+sig_type_call])
	message:
		"Drawing distribution plots for genome... [{wildcards.sampleName}]"

	shell:
		"""
		module load python3/3.6.3
		rBis.draw_dist_plots.py \
		-c {input.cov} -o {sample_dir}/{wildcards.sampleName}/{genome_dir}/call_{coverage_filtered_dir} \
		-n {wildcards.sampleName} \
		-s {sig_thresh_call} \
		-t {sig_type_call}
		"""

rule align_se_circRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + genome_dir + "/Align/unaligned.fq.gz"
	output:
		summary = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.log",
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.bam.bai",
		unaligned = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/unaligned_pre.fq.gz"
	message:
		"Aligning circRNA... [{wildcards.sampleName}]"
	threads:
		cluster["align_se"]["cpu"]
	params:
		desDir = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align"
	shell:
		"""
		module purge
		module load hisat2/2.2.1
		module load samtools/1.18.0

		mkdir -p {params.desDir}

		hisat-3n --un-gz {output.unaligned} {soft_clip} {RNA_splice} --no-unal --summary-file {output.summary} \
		-x {circRNA_index} \
		-q -U {input.fq1} \
		--base-change C,T --repeat --directional-mapping | samtools sort -o {output.bam}

		samtools index {output.bam}
		"""

rule get_align_stat_circRNA:
	input:
		summary = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.log"
	output:
		alignStats = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/alignStats.tsv"
	message:
		"Getting alignment stats... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.alignStat_perSample.py -s {input.summary} -o {output.alignStats} -n {wildcards.sampleName}
		"""

rule post_process_bam_circRNA:
	input:
		fq1 = sample_dir + "/{sampleName}/" + genome_dir + "/Align/unaligned.fq.gz",
		unalignedFq = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/unaligned_pre.fq.gz",
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.bam.bai",
	output:
		newUnaligned = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/unaligned.fq.gz",
		plusBam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam.bai",
		plusUniqBam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam",
		plusUniqBamBai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.minus.bam.bai",
	message:
		"Post-processing bam file from circRNA alignment... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0
		module load seqtk/1.3

		rBis.postProcessAlign.sh \
		-o {sample_dir}/{wildcards.sampleName}/{circRNA_dir}/Align \
		{input.bam} {input.fq1} {input.unalignedFq}

		"""

rule feature_count_all_circRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + circRNA_dir + "/featureCount_all.txt"
	message:
		"Feature count all circRNA_dir... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.count_chromosomes.py -b {input.bam} > {output.cnt}
		"""

rule feature_count_uniq_circRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam.bai",
	output:
		cnt = sample_dir + "/{sampleName}/" + circRNA_dir + "/featureCount_uniq.txt"
	message:
		"Feature count uniq piRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3
		rBis.count_chromosomes.py -b {input.bam} > {output.cnt}
		"""

#get candidates files for sample
def get_all_featureCounts_for_sample(sampleName):
	
	#get candidate files for sample and add them to list
	tRNA_all = sample_dir + "/" + sampleName + "/" + tRNA_dir + "/featureCount_all.txt"
	tRNA_uniq = sample_dir + "/" + sampleName + "/" + tRNA_dir + "/featureCount_uniq.txt"
	miRNA_all = sample_dir + "/" + sampleName + "/" + miRNA_dir + "/featureCount_all.txt"
	miRNA_uniq = sample_dir + "/" + sampleName + "/" + miRNA_dir + "/featureCount_uniq.txt"
	piRNA_all = sample_dir + "/" + sampleName + "/" + piRNA_dir + "/featureCount_all.txt"
	piRNA_uniq = sample_dir + "/" + sampleName + "/" + piRNA_dir + "/featureCount_uniq.txt"
	genome_all = sample_dir + "/" + sampleName + "/" + genome_dir + "/featureCount_all.txt"
	genome_uniq = sample_dir + "/" + sampleName + "/" + genome_dir + "/featureCount_uniq.txt"
	circRNA_all = sample_dir + "/" + sampleName + "/" + circRNA_dir + "/featureCount_all.txt"
	circRNA_uniq = sample_dir + "/" + sampleName + "/" + circRNA_dir + "/featureCount_uniq.txt"
	return [tRNA_all, tRNA_uniq, miRNA_all, miRNA_uniq, piRNA_all, piRNA_uniq, genome_all, genome_uniq, circRNA_all, circRNA_uniq]


rule merge_featureCounts_by_sample:
	input:
		featureCounts = lambda wildcards: get_all_featureCounts_for_sample(wildcards.sampleName)
	output:
		out = sample_dir + "/{sampleName}/featureCounts.tsv"
	message:
		"Merging all featureCounts for [{wildcards.sampleName}]..."
	shell:
		"""
		module load python3/3.6.3
		rBis.merge_featureCounts_by_sample.py \
		-f {input.featureCounts} \
		-o {sample_dir}/{wildcards.sampleName}/featureCounts.tsv \
		-s {wildcards.sampleName}
		"""

rule get_read_stats_circRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.bam",
		plusBam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam",
		plusBamBai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam.bai",
		minusBam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.minus.bam",
		minusBamBai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.minus.bam.bai",
		unalignedFQ = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/unaligned.fq.gz"

	output:
		readStats = sample_dir + "/{sampleName}/" + circRNA_dir + "/readStats.tsv"
	message:
		"Counting reads used in mpileup for circRNA... [{wildcards.sampleName}]"
	shell:
		"""
		module load samtools/1.18.0
		rBis.getReadStats.sh -o {output.readStats} {input.bam} {input.plusBam} {input.minusBam} {wildcards.sampleName} {input.unalignedFQ}
		"""

rule run_samtools_mpileup_all_circRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + circRNA_dir + "/mpileup_out_raw_all.tsv"
	message:
		"Running samtools mpileup for circRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {circRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}
		"""

rule run_samtools_mpileup_uniq_circRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam.bai",
	output:
		table = sample_dir + "/{sampleName}/" + circRNA_dir + "/mpileup_out_raw_uniq.tsv"
	message:
		"Running samtools mpileup uniq for circRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load samtools/1.18.0

		samtools mpileup -B -f {circRNA_fa} \
		--no-output-ins --no-output-del --no-output-ends --max-depth 10000 --min-BQ 30 --excl-flags 1540 \
		{input.bam} \
		| gawk '{{ if ( $3=="C" || $3=="c" || $3=="G" || $3=="g") print }}' \
		| gawk -F "\t" -v OFS="\t" '{{ gsub(/[<>]/, "", $5); gsub(/[-]/, "", $5); gsub(/[1-9]/, "", $5); print }}' \
		| gawk -F'\t' '{{ if ($5!="") print }}' \
		| gawk -F'\t' '{{ if ($4!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
		> {output.table}

		sed  -i '1i ref\tpos\trefBase\tcov\treadBase' {output.table}

		"""

rule call_m5C_candidates_circRNA:
	input:
		tableAll = sample_dir + "/{sampleName}/" + circRNA_dir + "/mpileup_out_raw_all.tsv",
		tableUniq = sample_dir + "/{sampleName}/" + circRNA_dir + "/mpileup_out_raw_uniq.tsv"
	output:
		cov = sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv",
		stats = sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv",
		statsSimplified = sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv"
	params:
		cov = cov_thresh,
		desDir = sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir
	message:
		"Calling m5C candidates for circRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module purge
		module load anaconda3
		source activate BisKit

		mkdir -p {params.desDir}

		rBis.mpileup_getStrand.py -c {input.tableAll} -o $TMPDIR/{wildcards.sampleName}_mpileup_all_circRNA.tsv
		rBis.mpileup_getStrand.py -c {input.tableUniq} -o $TMPDIR/{wildcards.sampleName}_mpileup_uniq_circRNA.tsv

		rBis.mpileup_mergeTables.py \
		-a $TMPDIR/{wildcards.sampleName}_mpileup_all_circRNA.tsv \
		-u $TMPDIR/{wildcards.sampleName}_mpileup_uniq_circRNA.tsv \
		-o $TMPDIR/{wildcards.sampleName}_mpileup_final_circRNA.tsv \

		rBis.mpileup_call_m5C.py \
		-v {params.cov} -c $TMPDIR/{wildcards.sampleName}_mpileup_final_circRNA.tsv \
		-o {sample_dir}/{wildcards.sampleName}/{circRNA_dir}/call_{coverage_filtered_dir} \
		-m circRNA -s {sig_thresh_call} -t {sig_type_call} -l {lookup_table} -n {min_C_count} \
		-e {wildcards.sampleName} -u {species} -r {genome_fa}

		rm $TMPDIR/{wildcards.sampleName}_mpileup_all_circRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_uniq_circRNA.tsv
		rm $TMPDIR/{wildcards.sampleName}_mpileup_final_circRNA.tsv

		"""

rule draw_distribution_plots_circRNA:
	input:
		cov = sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	output:
		plots = expand(sample_dir + "/{{sampleName}}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/{pre}.png", pre=["log2_cov", "methRate", "log2_cov_"+sig_type_call, "methRate_"+sig_type_call])
	message:
		"Drawing distribution plots for circRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module load python3/3.6.3

		rBis.draw_dist_plots.py \
		-c {input.cov} -o {sample_dir}/{wildcards.sampleName}/{circRNA_dir}/call_{coverage_filtered_dir} \
		-n {wildcards.sampleName} \
		-s {sig_thresh_call} \
		-t {sig_type_call}
		"""

rule merge_all_featureCounts:
	input:
		featureCounts = expand( sample_dir + "/{sampleName}/featureCounts.tsv" , sampleName=samples.Name.tolist())
	output:
		out = sample_dir + "/featureCounts.tsv"
	message:
		"Merging all featureCounts..."
	shell:
		"""
		module load python3/3.6.3
		rBis.merge_all_featureCounts.py \
		-f {input.featureCounts} \
		-o {sample_dir}/featureCounts.tsv
		"""


#get candidates files for sample
def get_all_alignStats_for_sample(sampleName):
	
	#get candidate files for sample and add them to list
	rRNA = sample_dir + "/" + sampleName + "/" + rRNA_dir + "/Align/alignStats.tsv" 
	tRNA = sample_dir + "/" + sampleName + "/" + tRNA_dir + "/Align/alignStats.tsv" 
	miRNA = sample_dir + "/" + sampleName + "/" + miRNA_dir + "/Align/alignStats.tsv"
	piRNA = sample_dir + "/" + sampleName + "/" + piRNA_dir + "/Align/alignStats.tsv"
	genome = sample_dir + "/" + sampleName + "/" + genome_dir + "/Align/alignStats.tsv"
	circRNA = sample_dir + "/" + sampleName + "/" + circRNA_dir + "/Align/alignStats.tsv" 
	return [rRNA, tRNA, miRNA, piRNA, genome, circRNA]

rule merge_align_stats:
	input:
		alignStats = lambda wildcards: get_all_alignStats_for_sample(wildcards.sampleName)
	output:
		mergedAlignStat = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedAlignStats.tsv"
	message:
		"Merging alignment statistics per sample... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3

		rBis.mergeAlignStats.py \
		-a {input.alignStats} \
		-o {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/mergedAlignStats \
		-s {wildcards.sampleName}
		"""

rule merge_align_stats_for_all_samples:
	input:
		rRNAstats = expand( sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/alignStats.tsv" , sampleName=samples.Name.tolist()),
		tRNAstats = expand( sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/alignStats.tsv" , sampleName=samples.Name.tolist()),
		miRNAstats = expand( sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/alignStats.tsv" , sampleName=samples.Name.tolist()),
		piRNAstats = expand( sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/alignStats.tsv" , sampleName=samples.Name.tolist()),
		genomeStats = expand( sample_dir + "/{sampleName}/" + genome_dir + "/Align/alignStats.tsv" , sampleName=samples.Name.tolist()),
		circRNAstats = expand( sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/alignStats.tsv" , sampleName=samples.Name.tolist())
	output:
		rRNA_out = mqc_dir + "/rRNA_alignStats.tsv",
		tRNA_out = mqc_dir + "/tRNA_alignStats.tsv",
		miRNA_out = mqc_dir + "/miRNA_alignStats.tsv",
		piRNA_out = mqc_dir + "/piRNA_alignStats.tsv",
		genome_out = mqc_dir + "/Genome_alignStats.tsv",
		circRNA_out = mqc_dir + "/circRNA_alignStats.tsv"
	message:
		"Combining all alignment statistics by source..."
	shell:
		"""
		echo | awk -vz="Sample" -va="Uniquely Aligned" -vc="Multi-Aligned" -ve="Unaligned" '{{ print z"\t"a"\t"c"\t"e }}' > $TMPDIR/tmprRNA.tsv
		cat $TMPDIR/tmprRNA.tsv {input.rRNAstats} | sed -n 'p;n' > {output.rRNA_out}
		rm $TMPDIR/tmprRNA.tsv

		echo | awk -vz="Sample" -va="Uniquely Aligned" -vc="Multi-Aligned" -ve="Unaligned" '{{ print z"\t"a"\t"c"\t"e }}' > $TMPDIR/tmptRNA.tsv
		cat $TMPDIR/tmptRNA.tsv {input.tRNAstats} | sed -n 'p;n' > {output.tRNA_out}
		rm $TMPDIR/tmptRNA.tsv

		echo | awk -vz="Sample" -va="Uniquely Aligned" -vc="Multi-Aligned" -ve="Unaligned" '{{ print z"\t"a"\t"c"\t"e }}' > $TMPDIR/tmpmiRNA.tsv
		cat $TMPDIR/tmpmiRNA.tsv {input.miRNAstats} | sed -n 'p;n' > {output.miRNA_out}
		rm $TMPDIR/tmpmiRNA.tsv

		echo | awk -vz="Sample" -va="Uniquely Aligned" -vc="Multi-Aligned" -ve="Unaligned" '{{ print z"\t"a"\t"c"\t"e }}' > $TMPDIR/tmppiRNA.tsv
		cat $TMPDIR/tmppiRNA.tsv {input.piRNAstats} | sed -n 'p;n' > {output.piRNA_out}
		rm $TMPDIR/tmppiRNA.tsv

		echo | awk -vz="Sample" -va="Uniquely Aligned" -vc="Multi-Aligned" -ve="Unaligned" '{{ print z"\t"a"\t"c"\t"e }}' > $TMPDIR/tmpgenome.tsv
		cat $TMPDIR/tmpgenome.tsv {input.genomeStats} | sed -n 'p;n' > {output.genome_out}
		rm $TMPDIR/tmpgenome.tsv

		echo | awk -vz="Sample" -va="Uniquely Aligned" -vc="Multi-Aligned" -ve="Unaligned" '{{ print z"\t"a"\t"c"\t"e }}' > $TMPDIR/tmpcircRNA.tsv
		cat $TMPDIR/tmpcircRNA.tsv {input.circRNAstats} | sed -n 'p;n' > {output.circRNA_out}
		rm $TMPDIR/tmpcircRNA.tsv

		"""

rule draw_read_stratification:
	input:
		mergedAlignStat = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedAlignStats.tsv",
		rRNAreadStat = sample_dir + "/{sampleName}/" + rRNA_dir + "/readStats.tsv",
		tRNAreadStat = sample_dir + "/{sampleName}/" + tRNA_dir + "/readStats.tsv",
		miRNAreadStat = sample_dir + "/{sampleName}/" + miRNA_dir + "/readStats.tsv",
		piRNAreadStat = sample_dir + "/{sampleName}/" + piRNA_dir + "/readStats.tsv",
		circRNAreadStat = sample_dir + "/{sampleName}/" + circRNA_dir + "/readStats.tsv"
	output:
		expand(sample_dir + "/{{sampleName}}/" + coverage_filtered_dir + "/" + per_sample_plot_dir + "/Read_Stratification_With{ext}", ext=["out_Unaligned.pdf", "out_Unaligned.png", "_Unaligned.pdf", "_Unaligned.png"])
	message:
		"Drawing read stratification per sample... [{wildcards.sampleName}]"
	params:
		desDir = sample_dir + "/{sampleName}/" + coverage_filtered_dir + per_sample_plot_dir
	shell:
		"""

		module load R/4.1.1

		mkdir -p {params.desDir}

		rBis.drawReadStratification.r \
		-a {input.mergedAlignStat} \
		-r {input.rRNAreadStat} \
		-t {input.tRNAreadStat} \
		-m {input.miRNAreadStat} \
		-p {input.piRNAreadStat} \
		-c {input.circRNAreadStat} \
		-o {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/{per_sample_plot_dir}/Read_Stratification
		"""

rule draw_m5C_stratification:
	input:
		mergedCandidates = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCandidates.tsv"
	output:
		expand(sample_dir + "/{{sampleName}}/" + coverage_filtered_dir + "/" + per_sample_plot_dir + "/Candidate_Stratification{ext}", sampleName=sample_list, ext=["_(All).pdf", "_(All).png", "_(" + sig_type_call + ").pdf", "_(" + sig_type_call + ").png", "_(" + sig_type_call + ")_(Annotated).pdf", "_(" + sig_type_call + ")_(Annotated).png"])
	message:
		"Drawing m5C stratification per sample... [{wildcards.sampleName}]"
	params:
		desDir = sample_dir + "/{sampleName}/" + coverage_filtered_dir + per_sample_plot_dir
	shell:
		"""
		module load R/4.1.1

		mkdir -p {params.desDir}

		rBis.drawM5CstratBarPlotPerSample.r \
		-s {sig_thresh_call} -t {sig_type_call} -r {mr_thresh} -c {cov_thresh} \
		-o {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/{per_sample_plot_dir}/Candidate_Stratification {input.mergedCandidates}
		"""

rule draw_read_stratification_bar:
	input:
		alignStats = expand(sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedAlignStats.tsv", sampleName=samples.Name.tolist()),
		rRNAstats = expand(sample_dir + "/{sampleName}/" + rRNA_dir + "/readStats.tsv", sampleName=samples.Name.tolist()),
		tRNAstats = expand(sample_dir + "/{sampleName}/" + tRNA_dir + "/readStats.tsv", sampleName=samples.Name.tolist()),
		miRNAstats = expand(sample_dir + "/{sampleName}/" + miRNA_dir + "/readStats.tsv", sampleName=samples.Name.tolist()),
		piRNAstats = expand(sample_dir + "/{sampleName}/" + piRNA_dir + "/readStats.tsv", sampleName=samples.Name.tolist()),
		genomeStats = expand(sample_dir + "/{sampleName}/" + genome_dir + "/readStats.tsv", sampleName=samples.Name.tolist()),
		circRNAstats = expand(sample_dir + "/{sampleName}/" + circRNA_dir + "/readStats.tsv", sampleName=samples.Name.tolist())
	output:
		stack = expand(plot_dir + "/Read_Stratification_Count{ext}", ext=[".pdf", ".png"]),
		percent = expand(plot_dir + "/Read_Stratification_Percentage{ext}", ext=[".pdf", ".png"]),
		tables = expand(mqc_dir + "/Read_Stratification_{suf}.tsv", suf = ["All", "Aligned"])
	message:
		"Drawing read stratification for all samples..."
	shell:
		"""
		module load R/4.1.1

		mkdir -p {plot_dir} {mqc_dir}

		rBis.drawReadStratBarPlot.r \
		-o {plot_dir}/Read_Stratification \
		-r {mqc_dir}/Read_Stratification \
		{input.alignStats} {input.rRNAstats} {input.tRNAstats} {input.miRNAstats} {input.piRNAstats} {input.genomeStats} {input.circRNAstats}
		"""


rule draw_m5C_stratification_bar:
	input:
		mergedM5C = expand(sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCandidates.tsv", sampleName=samples.Name.tolist())
	output:
		allStack = expand(plot_dir + "/Candidate_Stratification_{ext}", ext=["Count_(All).pdf", "Count_(All).png", "Percentage_(All).pdf", "Percentage_(All).png"]),
		significant = expand(plot_dir + "/Candidate_Stratification_{ext}", ext=["Count_(Significant).pdf", "Count_(Significant).png", "Percentage_(Significant).pdf", "Percentage_(Significant).png"]),
		genome = expand(plot_dir + "/Candidate_Stratification_{ext}", ext=["Count_(Annotated).pdf", "Count_(Annotated).png", "Percentage_(Annotated).pdf", "Percentage_(Annotated).png"])
	message:
		"Drawing m5C stratification for all samples..."
	shell:
		"""
		module load R/4.1.1

		rBis.drawM5CstratBarPlotAllSamples.r \
		-s {sig_thresh_call} -t {sig_type_call} -r {mr_thresh} -c {cov_thresh} \
		-o {plot_dir}/Candidate_Stratification {input.mergedM5C}
		"""

#get candidates files for sample
def get_all_callStats_for_sample(sampleName):
	
	#get candidate files for sample and add them to list
	rRNA = sample_dir + "/" + sampleName + "/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv" 
	tRNA = sample_dir + "/" + sampleName + "/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv"
	miRNA = sample_dir + "/" + sampleName + "/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv"
	piRNA = sample_dir + "/" + sampleName + "/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv"
	genome = sample_dir + "/" + sampleName + "/" + genome_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv"
	circRNA = sample_dir + "/" + sampleName + "/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv"
	return [rRNA, tRNA, miRNA, piRNA, genome, circRNA]

rule merge_call_stats:
	input:
		callStats = lambda wildcards: get_all_callStats_for_sample(wildcards.sampleName)
	output:
		mergedCallStat = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCallStats.tsv"
	message:
		"Merging call statistics per sample... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3

		rBis.mergeCallStats.py \
		-a {input.callStats} \
		-o {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/mergedCallStats \
		-s {wildcards.sampleName}
		"""

rule merge_call_stats_for_all_samples:
	input:
		rRNAstats = expand( sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv", sampleName=samples.Name.tolist()),
		tRNAstats = expand( sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv", sampleName=samples.Name.tolist()),
		miRNAstats = expand( sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv", sampleName=samples.Name.tolist()),
		piRNAstats = expand( sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv", sampleName=samples.Name.tolist()),
		genomeStats = expand( sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv", sampleName=samples.Name.tolist()),
		circRNAstats = expand( sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats.tsv", sampleName=samples.Name.tolist()),
		rRNAstatsS = expand( sample_dir + "/{sampleName}/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv", sampleName=samples.Name.tolist()),
		tRNAstatsS = expand( sample_dir + "/{sampleName}/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv", sampleName=samples.Name.tolist()),
		miRNAstatsS = expand( sample_dir + "/{sampleName}/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv", sampleName=samples.Name.tolist()),
		piRNAstatsS = expand( sample_dir + "/{sampleName}/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv", sampleName=samples.Name.tolist()),
		genomeStatsS = expand( sample_dir + "/{sampleName}/" + genome_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv", sampleName=samples.Name.tolist()),
		circRNAstatsS = expand( sample_dir + "/{sampleName}/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call_stats_simplified.tsv", sampleName=samples.Name.tolist())
	output:
		rRNA_out = mqc_dir + "/rRNA_callStats.tsv",
		tRNA_out = mqc_dir + "/tRNA_callStats.tsv",
		miRNA_out = mqc_dir + "/miRNA_callStats.tsv",
		piRNA_out = mqc_dir + "/piRNA_callStats.tsv",
		genome_out = mqc_dir + "/Genome_callStats.tsv",
		circRNA_out = mqc_dir + "/circRNA_callStats.tsv",
		rRNA_updown = mqc_dir + "/rRNA_callStats_sigNonsig.tsv",
		tRNA_updown = mqc_dir + "/tRNA_callStats_sigNonsig.tsv",
		miRNA_updown = mqc_dir + "/miRNA_callStats_sigNonsig.tsv",
		piRNA_updown = mqc_dir + "/piRNA_callStats_sigNonsig.tsv",
		genome_updown = mqc_dir + "/Genome_callStats_sigNonsig.tsv",
		circRNA_updown = mqc_dir + "/circRNA_callStats_sigNonsig.tsv"
	message:
		"Combining all call statistics..."
	shell:
		"""
		mkdir -p {mqc_dir}
		
		echo | awk -vz="Sample" -vt="Total Candidates" -va=">= 10 cov" -vb=">= 0.1 MethRate" -vc="< {sig_thresh_call} {sig_type_call}" -vd="Significant Candidates" '{{ print z"\t"t"\t"a"\t"b"\t"c"\t"d }}' > $TMPDIR/$$_rRNACallStats.tsv
		cat $TMPDIR/$$_rRNACallStats.tsv {input.rRNAstats} | sed -n 'p;n' > {output.rRNA_out}

		echo | awk -vz="Sample" -va="Significant Candidates" -vb="Non-Significant Candidates" '{{ print z"\t"a"\t"b }}' > $TMPDIR/$$_rRNACallStatsS.tsv
		cat $TMPDIR/$$_rRNACallStatsS.tsv {input.rRNAstatsS} | sed -n 'p;n' > {output.rRNA_updown}
		rm $TMPDIR/$$_rRNACallStats.tsv
		rm $TMPDIR/$$_rRNACallStatsS.tsv

		echo | awk -vz="Sample" -vt="Total Candidates" -va=">= 10 cov" -vb=">= 0.1 MethRate" -vc="< {sig_thresh_call} {sig_type_call}" -vd="Significant Candidates" '{{ print z"\t"t"\t"a"\t"b"\t"c"\t"d }}' > $TMPDIR/$$_tRNACallStats.tsv
		cat $TMPDIR/$$_tRNACallStats.tsv {input.tRNAstats} | sed -n 'p;n' > {output.tRNA_out}

		echo | awk -vz="Sample" -va="Significant Candidates" -vb="Non-Significant Candidates" '{{ print z"\t"a"\t"b }}' > $TMPDIR/$$_tRNACallStatsS.tsv
		cat $TMPDIR/$$_tRNACallStatsS.tsv {input.tRNAstatsS} | sed -n 'p;n' > {output.tRNA_updown}
		rm $TMPDIR/$$_tRNACallStats.tsv
		rm $TMPDIR/$$_tRNACallStatsS.tsv

		echo | awk -vz="Sample" -vt="Total Candidates" -va=">= 10 cov" -vb=">= 0.1 MethRate" -vc="< {sig_thresh_call} {sig_type_call}" -vd="Significant Candidates" '{{ print z"\t"t"\t"a"\t"b"\t"c"\t"d }}' > $TMPDIR/$$_miRNACallStats.tsv
		cat $TMPDIR/$$_miRNACallStats.tsv {input.miRNAstats} | sed -n 'p;n' > {output.miRNA_out}

		echo | awk -vz="Sample" -va="Significant Candidates" -vb="Non-Significant Candidates" '{{ print z"\t"a"\t"b }}' > $TMPDIR/$$_miRNACallStatsS.tsv
		cat $TMPDIR/$$_miRNACallStatsS.tsv {input.miRNAstatsS} | sed -n 'p;n' > {output.miRNA_updown}
		rm $TMPDIR/$$_miRNACallStats.tsv
		rm $TMPDIR/$$_miRNACallStatsS.tsv

		echo | awk -vz="Sample" -vt="Total Candidates" -va=">= 10 cov" -vb=">= 0.1 MethRate" -vc="< {sig_thresh_call} {sig_type_call}" -vd="Significant Candidates" '{{ print z"\t"t"\t"a"\t"b"\t"c"\t"d }}' > $TMPDIR/$$_piRNACallStats.tsv
		cat $TMPDIR/$$_piRNACallStats.tsv {input.piRNAstats} | sed -n 'p;n' > {output.piRNA_out}

		echo | awk -vz="Sample" -va="Significant Candidates" -vb="Non-Significant Candidates" '{{ print z"\t"a"\t"b }}' > $TMPDIR/$$_piRNACallStatsS.tsv
		cat $TMPDIR/$$_piRNACallStatsS.tsv {input.piRNAstatsS} | sed -n 'p;n' > {output.piRNA_updown}
		rm $TMPDIR/$$_piRNACallStats.tsv
		rm $TMPDIR/$$_piRNACallStatsS.tsv

		echo | awk -vz="Sample" -vt="Total Candidates" -va=">= 10 cov" -vb=">= 0.1 MethRate" -vc="< {sig_thresh_call} {sig_type_call}" -vd="Significant Candidates" '{{ print z"\t"t"\t"a"\t"b"\t"c"\t"d }}' > $TMPDIR/$$_genomeCallStats.tsv
		cat $TMPDIR/$$_genomeCallStats.tsv {input.genomeStats} | sed -n 'p;n' > {output.genome_out}

		echo | awk -vz="Sample" -va="Significant Candidates" -vb="Non-Significant Candidates" '{{ print z"\t"a"\t"b }}' > $TMPDIR/$$_genomeCallStatsS.tsv
		cat $TMPDIR/$$_genomeCallStatsS.tsv {input.genomeStatsS} | sed -n 'p;n' > {output.genome_updown}
		rm $TMPDIR/$$_genomeCallStats.tsv
		rm $TMPDIR/$$_genomeCallStatsS.tsv

		echo | awk -vz="Sample" -vt="Total Candidates" -va=">= 10 cov" -vb=">= 0.1 MethRate" -vc="< {sig_thresh_call} {sig_type_call}" -vd="Significant Candidates" '{{ print z"\t"t"\t"a"\t"b"\t"c"\t"d }}' > $TMPDIR/$$_circRNACallStats.tsv
		cat $TMPDIR/$$_circRNACallStats.tsv {input.circRNAstats} | sed -n 'p;n' > {output.circRNA_out}

		echo | awk -vz="Sample" -va="Significant Candidates" -vb="Non-Significant Candidates" '{{ print z"\t"a"\t"b }}' > $TMPDIR/$$_circRNACallStatsS.tsv
		cat $TMPDIR/$$_circRNACallStatsS.tsv {input.circRNAstatsS} | sed -n 'p;n' > {output.circRNA_updown}
		rm $TMPDIR/$$_circRNACallStats.tsv
		rm $TMPDIR/$$_circRNACallStatsS.tsv

		"""

rule merge_all_call_stats_for_multiqc:
	input:
		allCandidates = expand(mqc_dir + "/{src}_callStats.tsv", src = ["rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA"]),
		sigNonsig = expand(mqc_dir + "/{src}_callStats_sigNonsig.tsv", src = ["rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA"])
	output:
		allCandidates = mqc_dir + "/callStats_all.tsv",
		sigNonsig = mqc_dir + "/callStats_sigNonsig.tsv"
	message:
		"Combining call stats..."
	shell:
		"""
		module load python3/3.6.3

		mkdir -p {mqc_dir}

		rBis.merge_all_callStats.py \
		-a {input.allCandidates} \
		-s {input.sigNonsig} \
		-o {mqc_dir}

		"""

#get candidates files for sample
def get_all_candidates_for_sample(sampleName):
	
	#get candidate files for sample and add them to list
	rRNA = sample_dir + "/" + sampleName + "/" + rRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	tRNA = sample_dir + "/" + sampleName + "/" + tRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	miRNA = sample_dir + "/" + sampleName + "/" + miRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	piRNA = sample_dir + "/" + sampleName + "/" + piRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	genome = sample_dir + "/" + sampleName + "/" + genome_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	circRNA = sample_dir + "/" + sampleName + "/" + circRNA_dir + "/call_" + coverage_filtered_dir + "/call.tsv"
	return [rRNA, tRNA, miRNA, piRNA, genome, circRNA]


## merge m5C candidates by each sample
rule merge_call_results:
	input:
		candidates = lambda wildcards: get_all_candidates_for_sample(wildcards.sampleName)
	output:
		mergedCandidates = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCandidates.tsv"
	message:
		"Merging m5C candidates per sample... [{wildcards.sampleName}]"
	shell:
		"""
		module load python3/3.6.3



		rBis.mergeCandidates.py \
		-c {input.candidates} \
		-o {output.mergedCandidates}
		"""

## merge m5C candidates by each sample
rule draw_candidate_profile:
	input:
		mergedCandidates = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCandidates.tsv"
	output:
		plusBW = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/plus.bw",
		minusBW = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/minus.bw"
	message:
		"drawing candidate profile... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery

		rBis.draw_candidate_profile.py \
		-c {input.mergedCandidates} \
		-o {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir} \
		-l {lookup_table}

		bedGraphToBigWig {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/plusCandidates.bg {chrom_size} {output.plusBW}
		bedGraphToBigWig {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/minusCandidates.bg {chrom_size} {output.minusBW}

		rm {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/plusCandidates.bg
		rm {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/minusCandidates.bg
		"""

#get candidates files for sample
def get_comparison_candidate_tsv(diffPairName):

	#get candidate files for sample and add them to list
	groups = config_pw['diff_pair_list'][diffPairName]
	currentSamples = []
	for group in groups:
		samp = samples.Name[samples.Group == group].tolist()
		currentSamples += samp

	mergedCandidateList = []

	for samp in currentSamples:
		mergedCandidateTSV = sample_dir + "/" + samp + "/" + coverage_filtered_dir + "/mergedCandidates.tsv"
		mergedCandidateList.append(mergedCandidateTSV)

	return mergedCandidateList


## compare all candidates per comparison
rule compare_candidates:
	input:
		candidates = lambda wildcards: get_comparison_candidate_tsv(wildcards.diffPairName)
	output:
		tsv = compare_dir + "/{diffPairName}/pairwise_comparison.tsv"
	params:
		groupList = lambda wildcards: get_group(wildcards.diffPairName),
		desDir = compare_dir + "/{diffPairName}"
	message:
		"Comparing candidates between two samples... [{wildcards.diffPairName}]"
	shell:
		"""
		module load python3/3.6.3

		mkdir -p {params.desDir}

		rBis.compare_pairwise.py -c {input.candidates} -o {compare_dir}/{wildcards.diffPairName} -s {src_sampleInfo} -p {wildcards.diffPairName} -g {params.groupList} -m {min_reps}
		"""

def get_group(diffPairName):
	groupList = diffPairNameDict[diffPairName]
	return ','.join(groupList)

## categorize by threshold
rule categorize_comparisons:
	input:
		pairwise = compare_dir + "/{diffPairName}/pairwise_comparison.tsv"
	output:
		tsv = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/categorized_pairwise_comparison.tsv",
		stats = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/stats.tsv",
		plots = expand(compare_dir + "/{{diffPairName}}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/{pre}.png", pre=["methRate_log2foldChange", "methRate_log2foldChange_"+sig_type_diff, "delta_methRate", "delta_methRate_"+sig_type_diff])
	params:
		sig = sig_thresh_diff,
		diff = diff_thresh,
		groupList = lambda wildcards: get_group(wildcards.diffPairName),
		desDir = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh
	message:
		"Categorizing comparison between two samples... [{wildcards.diffPairName}]"
	shell:
		"""
		module load python3/3.6.3

		mkdir -p {params.desDir}

		rBis.categorize_comparisons.py \
		-p {input.pairwise} \
		-o {compare_dir}/{wildcards.diffPairName}/Cov{cov_thresh}_{sig_type_diff}{sig_thresh_diff}_Diff{diff_thresh} \
		-s {sig_thresh_diff} -t {sig_type_diff} \
		-d {params.diff} \
		-i {src_sampleInfo} \
		-n {wildcards.diffPairName} \
		-g {params.groupList}
		"""


## get final comparisons by merging all pairwise sample comparisons
rule merge_all_pairwise:
	input:
		pairwises = expand(compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/categorized_pairwise_comparison.tsv", diffPairName=diffPairNameL),
	output:
		combined_compare_dir + "/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/all_sample_comparisons.tsv"
	message:
		"Combining all sample comparisons..."
	shell:
		"""
		module load python3/3.6.3
		rBis.merge_all_pairwise.py \
		-p {input.pairwises} \
		-o {combined_compare_dir}/Cov{cov_thresh}_{sig_type_diff}{sig_thresh_diff}_Diff{diff_thresh}
		"""


rule draw_sig_m5C_proportion_per_sample:
	input:
		mergedCallStats = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCallStats.tsv"
	output:
		expand(sample_dir + "/{{sampleName}}/" + coverage_filtered_dir + "/" + per_sample_plot_dir + "/Significant_Candidate_{ext}", sampleName=sample_list, ext=["Percentage.pdf", "Percentage.png", "Count.pdf", "Count.png"])
	message:
		"Drawing significant vs all m5C stats... [{wildcards.sampleName}]"
	params:
		desDir = sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/" + per_sample_plot_dir
	shell:
		"""
		module load R/4.1.1

		mkdir -p {params.desDir}

		rBis.drawM5CsigBarPlotPerSample.r \
		-o {sample_dir}/{wildcards.sampleName}/{coverage_filtered_dir}/{per_sample_plot_dir}/Significant_Candidate {input.mergedCallStats} \
		-t {sig_type_call}
		"""

rule draw_sig_m5C_proportion_all_samples:
	input:
		mergedCallStats = expand(sample_dir + "/{sampleName}/" + coverage_filtered_dir + "/mergedCallStats.tsv", sampleName=samples.Name.tolist()),
	output:
		sig = expand(plot_dir + "/Significant_Candidate_{ext}", ext=["Count.pdf", "Count.png", "Percentage.pdf", "Percentage.png"]),
	message:
		"Drawing significant vs all m5C stats for all samples..."
	shell:
		"""
		module load R/4.1.1

		mkdir -p {plot_dir}

		rBis.drawM5CsigBarPlotAllSample.r \
		-o {plot_dir}/Significant_Candidate {input.mergedCallStats} -t {sig_type_call}
		"""

rule draw_categorizations_by_source_per_pairwise:
	input:
		categorization = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/categorized_pairwise_comparison.tsv"
	output:
		plots = expand(compare_dir + "/{{diffPairName}}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/Categorization_By_Source_{ext}_{cat}", diffPairName=diffPairNameL,  ext=["Count", "Percentage"], cat=["(UPvDOWN).png", "(UPvDOWN).pdf", "(UNIQ1v2).png", "(UNIQ1v2).pdf", "(All).png", "(All).pdf", "(Single).png", "(Single).pdf"]),
		table = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/Categorization_By_Source.tsv"
	message:
		"Categorizing comparison between two samples... [{wildcards.diffPairName}]"
	params:
		desDir = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir
	shell:
		"""
		module load R/4.1.1
		rBis.drawCategorizationBySourcePerPairwise.r \
		-o {compare_dir}/{wildcards.diffPairName}/Cov{cov_thresh}_{sig_type_diff}{sig_thresh_diff}_Diff{diff_thresh}/{per_pairwise_plot_dir}/Categorization_By_Source {input.categorization}
		"""

rule combine_all_categorization_tables:
	input:
		categorization = expand(compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/Categorization_By_Source.tsv", diffPairName=diffPairNameL)
	output:
		table = combined_compare_dir + "/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/Categorization_By_Source_For_All_Pairwise_Comparisons.tsv"
	message:
		"Combining all pairwise categorizations..."
	shell:
		"""
		echo | awk -vs="Comparisons" -va="UP" -vb="DOWN" -ve="UNCHANGED" -vc="uniq1" -vd="uniq2" '{{ print s"\t"a"\t"b"\t"e"\t"c"\t"d }}' > $TMPDIR/tmp.tsv
		cat $TMPDIR/tmp.tsv {input.categorization} > {output.table}
		rm $TMPDIR/tmp.tsv
		"""

rule draw_categorizations_pairwise:
	input:
		categorizationStats = expand(compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/stats.tsv", diffPairName=diffPairNameL)
	output:
		plots = expand(combined_compare_dir + "/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/Differential_Analysis_Categorization_{ext}_{cat}", ext=["Count", "Percentage"], cat=["(UPvDOWN).png", "(UPvDOWN).pdf", "(UNIQ1v2).png", "(UNIQ1v2).pdf", "(All).png", "(All).pdf"]),
		table = mqc_dir + "/Differential_Analysis_Categorization.tsv"
	message:
		"Categorizing comparison between two samples..."
	shell:
		"""
		module load R/4.1.1
		rBis.drawCategorizationStatsPairwise.r \
		-o {combined_compare_dir}/Cov{cov_thresh}_{sig_type_diff}{sig_thresh_diff}_Diff{diff_thresh}/Differential_Analysis_Categorization -m {mqc_dir}/Differential_Analysis_Categorization {input.categorizationStats}
		"""

rule draw_categorizations_per_pairwise_as_volcano:
	input:
		categorization = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/categorized_pairwise_comparison.tsv"
	output:
		plots = expand(compare_dir + "/{{diffPairName}}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/volcano{ext}", ext=[".png", ".pdf"])
	message:
		"Categorizing comparison between two samples... [{wildcards.diffPairName}]"
	params:
		desDir = compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir
	shell:
		"""
		module load python3/3.6.3

		mkdir -p {params.desDir}

		rBis.drawVolcanoPerComparison.py \
		-c {input.categorization} \
		-t {sig_type_diff} \
		-o {compare_dir}/{wildcards.diffPairName}/Cov{cov_thresh}_{sig_type_diff}{sig_thresh_diff}_Diff{diff_thresh}/{per_pairwise_plot_dir}/volcano \
		-y 30
		"""

rule merge_volcano_for_mqc:
	input:
		expand(compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/volcano.png", diffPairName=diffPairNameL)
	output:
		mqc_dir + "/Volcano_Plots.png"
	message:
		"merging volcano plots for MultiQC Report visualization..."
	shell:
		"""
		module purge
		module load ImageMagick/7.0.7

		mkdir -p {mqc_dir}

		convert {input} -append {output}
		"""

rule separate_and_merge_categorizations_by_source:
	input:
		categorization = expand(compare_dir + "/{diffPairName}/Cov"+ cov_thresh + "_" + sig_type_diff + sig_thresh_diff + "_Diff" + diff_thresh + "/" + per_pairwise_plot_dir + "/Categorization_By_Source.tsv", diffPairName=diffPairNameL)
	output:
		outFiles = expand(mqc_dir + "/{src}_catStats.tsv", src=["rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA"]),
		outFilesUD = expand(mqc_dir + "/{src}_catStats_updown.tsv", src=["rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA"])
	message:
		"re-merging categorization stats..."
	shell:
		"""
		module load python3/3.6.3

		mkdir -p {mqc_dir}

		rBis.mergeCatStatsBySource.py \
		-p {input.categorization} \
		-o {mqc_dir}
		"""

rule get_deltaMR_dist_lines:
	input:
		pairwise = expand(compare_dir + "/{diffPairName}/pairwise_comparison.tsv", diffPairName=diffPairNameL)
	output:
		coords = expand(mqc_dir + "/delta_MR_{{diffPairName}}.tsv")
	message:
		"Create coordinates for delta methRate distribution line plot for MultiQC..."
	shell:
		"""
		module load R/4.1.1

		mkdir -p {mqc_dir}

		rBis.get_deltaMR_coords_for_lineplot.r \
		-o {mqc_dir} {input.pairwise}
		"""

rule create_multiqc:
	input:
		alignStats = expand(mqc_dir + "/{src}_alignStats.tsv", src= ["rRNA", "tRNA", "miRNA", "piRNA", "circRNA", "Genome"]),
		callStats = expand(mqc_dir + "/{src}_callStats.tsv", src= ["rRNA", "tRNA", "miRNA", "piRNA", "circRNA", "Genome"]),
		catStats = expand(mqc_dir + "/{src}_catStats.tsv", src= ["rRNA", "tRNA", "miRNA", "piRNA", "circRNA", "Genome"]),
		allCallStats = mqc_dir + "/callStats_all.tsv",
		allCallStatsSig = mqc_dir + "/callStats_sigNonsig.tsv",
		readStratification = expand(mqc_dir + "/Read_Stratification_{suf}.tsv", suf=["All", "Aligned"]),
		DifferentialAnalysisCategorization = mqc_dir + "/Differential_Analysis_Categorization.tsv",
		delta_MR_coords = expand(mqc_dir + "/delta_MR_{diffPairName}.tsv", diffPairName=diffPairNameL),
		m_bias = expand(mqc_dir + "/m-bias_{src}_{sampleName}.tsv", src= ["rRNA", "tRNA", "miRNA", "piRNA", "circRNA", "genome"], sampleName=samples.Name.tolist()),

	output:
		"multiqc_report.html"

	message:
		"Create MultiQC Report..."

	singularity:
		"docker://ewels/multiqc:v1.14"

	shell:
		"""
		multiqc {mqc_dir} -f -c {mqc_yaml}

		"""

rule get_readLen_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai"

	output:
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"
	message:
		"Creating m-bias plot for genome... [{wildcards.sampleName}]"
	params:
		desDir = sample_dir + "/{sampleName}/" + genome_dir
	shell:
		"""
		module purge
		module load samtools/1.18.0

		mkdir -p {params.desDir}

		## get most common read length
		#samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 > {output.readLenGenome}

		## get longest read length
		samtools stats {input.bam} | grep ^RL | cut -f 2- | tail -n 1 | cut -f1 > {output.readLenGenome}

		"""

rule draw_mbias_plot_genome:
	input:
		bam = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam",
		bai = sample_dir + "/{sampleName}/" + genome_dir + "/Align/align.bam.bai",
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"

	output:
		plots_commonReadLen = expand(sample_dir + "/{{sampleName}}/" + genome_dir + "/Align/m-bias_mostCommonReadLength_{direction}.png", direction = ["forward", "reverse"]),
		plots_commonReadLenGenome = expand(sample_dir + "/{{sampleName}}/" + genome_dir + "/Align/m-bias_mostCommonReadLengthOfGenomeAlignment_{direction}.png", direction = ["forward", "reverse"]),
		mqc_mbias = mqc_dir + "/m-bias_genome_{sampleName}.tsv"
	message:
		"Creating m-bias plot for genome... [{wildcards.sampleName}]"

	shell:
		"""
		module load python3/3.6.3
		module load samtools/1.18.0

		mkdir -p {mqc_dir}

		## get most common read length
		readLenCommon=$( samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 )

		readLenGenome=$( cat {input.readLenGenome} )

		rBis.draw_mBias_plot.py -b {input.bam} -r $readLenCommon -g $readLenGenome -o {sample_dir}/{wildcards.sampleName}/{genome_dir}/Align/m-bias > {output.mqc_mbias}

		"""

rule draw_mbias_plot_rRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + rRNA_dir + "/Align/align.plus.uniq.bam.bai",
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"

	output:
		plots_commonReadLen = expand(sample_dir + "/{{sampleName}}/" + rRNA_dir + "/Align/m-bias_mostCommonReadLength_{direction}.png", direction = ["forward", "reverse"]),
		plots_commonReadLenGenome = expand(sample_dir + "/{{sampleName}}/" + rRNA_dir + "/Align/m-bias_mostCommonReadLengthOfGenomeAlignment_{direction}.png", direction = ["forward", "reverse"]),
		mqc_mbias = mqc_dir + "/m-bias_rRNA_{sampleName}.tsv"
	message:
		"Creating m-bias plot for rRNA... [{wildcards.sampleName}]"

	shell:
		"""

		module load python3/3.6.3
		module load samtools/1.18.0

		mkdir -p {mqc_dir}

		## get most common read length
		readLenCommon=$( samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 )

		readLenGenome=$( cat {input.readLenGenome} )

		rBis.draw_mBias_plot.py -b {input.bam} -r $readLenCommon -g $readLenGenome -o {sample_dir}/{wildcards.sampleName}/{rRNA_dir}/Align/m-bias > {output.mqc_mbias}

		"""

rule draw_mbias_plot_tRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + tRNA_dir + "/Align/align.plus.uniq.bam.bai",
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"

	output:
		plots_commonReadLen = expand(sample_dir + "/{{sampleName}}/" + tRNA_dir + "/Align/m-bias_mostCommonReadLength_{direction}.png", direction = ["forward", "reverse"]),
		plots_commonReadLenGenome = expand(sample_dir + "/{{sampleName}}/" + tRNA_dir + "/Align/m-bias_mostCommonReadLengthOfGenomeAlignment_{direction}.png", direction = ["forward", "reverse"]),
		mqc_mbias = mqc_dir + "/m-bias_tRNA_{sampleName}.tsv"
	message:
		"Creating m-bias plot for tRNA... [{wildcards.sampleName}]"
	shell:
		"""

		module load python3/3.6.3
		module load samtools/1.18.0

		mkdir -p {mqc_dir}

		## get most common read length
		readLenCommon=$( samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 )

		readLenGenome=$( cat {input.readLenGenome} )

		rBis.draw_mBias_plot.py -b {input.bam} -r $readLenCommon -g $readLenGenome -o {sample_dir}/{wildcards.sampleName}/{tRNA_dir}/Align/m-bias > {output.mqc_mbias}

		"""

rule draw_mbias_plot_miRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + miRNA_dir + "/Align/align.plus.uniq.bam.bai",
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"

	output:
		plots_commonReadLen = expand(sample_dir + "/{{sampleName}}/" + miRNA_dir + "/Align/m-bias_mostCommonReadLength_{direction}.png", direction = ["forward", "reverse"]),
		plots_commonReadLenGenome = expand(sample_dir + "/{{sampleName}}/" + miRNA_dir + "/Align/m-bias_mostCommonReadLengthOfGenomeAlignment_{direction}.png", direction = ["forward", "reverse"]),
		mqc_mbias = mqc_dir + "/m-bias_miRNA_{sampleName}.tsv"
	message:
		"Creating m-bias plot for miRNA... [{wildcards.sampleName}]"

	shell:
		"""
		module load python3/3.6.3
		module load samtools/1.18.0

		mkdir -p {mqc_dir}

		## get most common read length
		readLenCommon=$( samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 )

		readLenGenome=$( cat {input.readLenGenome} )

		rBis.draw_mBias_plot.py -b {input.bam} -r $readLenCommon -g $readLenGenome -o {sample_dir}/{wildcards.sampleName}/{miRNA_dir}/Align/m-bias > {output.mqc_mbias}

		"""

rule draw_mbias_plot_piRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + piRNA_dir + "/Align/align.plus.uniq.bam.bai",
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"

	output:
		plots_commonReadLen = expand(sample_dir + "/{{sampleName}}/" + piRNA_dir + "/Align/m-bias_mostCommonReadLength_{direction}.png", direction = ["forward", "reverse"]),
		plots_commonReadLenGenome = expand(sample_dir + "/{{sampleName}}/" + piRNA_dir + "/Align/m-bias_mostCommonReadLengthOfGenomeAlignment_{direction}.png", direction = ["forward", "reverse"]),
		mqc_mbias = mqc_dir + "/m-bias_piRNA_{sampleName}.tsv"
	message:
		"Creating m-bias plot for piRNA... [{wildcards.sampleName}]"

	shell:
		"""

		module load python3/3.6.3
		module load samtools/1.18.0

		mkdir -p {mqc_dir}

		readLenGenome=$( cat {input.readLenGenome} )

		## get number of reads in uniq.plus.bam
		reads=$( samtools view -c {input.bam} )

		if [ $reads -eq 0 ]; then
			readLenCommon=$readLenGenome
		else
			## get most common read length
			readLenCommon=$( samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 )
		fi

		rBis.draw_mBias_plot.py -b {input.bam} -r $readLenCommon -g $readLenGenome -o {sample_dir}/{wildcards.sampleName}/{piRNA_dir}/Align/m-bias > {output.mqc_mbias}

		"""

rule draw_mbias_plot_circRNA:
	input:
		bam = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam",
		bai = sample_dir + "/{sampleName}/" + circRNA_dir + "/Align/align.plus.uniq.bam.bai",
		readLenGenome = sample_dir + "/{sampleName}/" + genome_dir + "/readLenGenome.txt"

	output:
		plots_commonReadLen = expand(sample_dir + "/{{sampleName}}/" + circRNA_dir + "/Align/m-bias_mostCommonReadLength_{direction}.png", direction = ["forward", "reverse"]),
		plots_commonReadLenGenome = expand(sample_dir + "/{{sampleName}}/" + circRNA_dir + "/Align/m-bias_mostCommonReadLengthOfGenomeAlignment_{direction}.png", direction = ["forward", "reverse"]),
		mqc_mbias = mqc_dir + "/m-bias_circRNA_{sampleName}.tsv"
	message:
		"Creating m-bias plot for circRNA... [{wildcards.sampleName}]"

	shell:
		"""

		module load python3/3.6.3
		module load samtools/1.18.0

		mkdir -p {mqc_dir}

		## get most common read length
		readLenCommon=$( samtools stats {input.bam} | grep ^RL | cut -f 2- | gawk 'BEGIN {{ max = 0; row = "" }} $2 > max {{ max = $2; row = $0 }} END {{ print row }}' | cut -f1 )

		readLenGenome=$( cat {input.readLenGenome} )

		rBis.draw_mBias_plot.py -b {input.bam} -r $readLenCommon -g $readLenGenome -o {sample_dir}/{wildcards.sampleName}/{circRNA_dir}/Align/m-bias > {output.mqc_mbias}

		"""