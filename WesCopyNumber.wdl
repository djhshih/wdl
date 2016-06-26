workflow WesCopyNumber {
	String tumor_name
	String normal_name
	File tumor_bam
	File normal_bam

	call runGatk4cnv { input:
		sample_name = tumor_name,
		sample_bam = tumor_bam
	}

	call callSnvsOnPair { input:
		tumor_name = tumor_name,
		normal_name = normal_name,
		tumor_bam = tumor_bam,
		normal_bam = normal_bam
	}

	call selectNormalHetSites { input:
		in = callSnvsOnPair.call_stats,
		out_fname = "${tumor_name}_germ-hets.call_stats"
	}	
		
	call convertCallStatsToCov { input:
		in = selectNormalHetSites.out,
		out_fname = "${tumor_name}_germ-hets.cov"
	}

	call runAllelicCapseg { input:
		sample_name = tumor_name,
		copy_ratio = runGatk4cnv.copy_ratio,
		seg = runGatk4cnv.seg,
		cov = convertCallStatsToCov.out
	}
}

task runGatk4cnv {
	String sample_name
	File sample_bam
	File sample_bai = sub(sample_bam, "\\.bam$", ".bai")
	File ref_fasta
	File ref_fai = ref_fasta + ".fai"
	File ref_dict = sub(ref_fasta, "\\.fasta$", ".dict")
	File  targets_bed
	Int padding
	File pon

	command {
		ln -s ${sample_bam} ${sample_name}.bam
		ln -s ${sample_bai} ${sample_name}.bai
		echo '${sample_name}.bam' > sample_bam.txt
		gatk4cnv \
			-i sample_bam.txt \
			--reference ${ref_fasta} \
			--intervalfile ${targets_bed} \
			--padding ${padding} \
			--keepdups \
			--ponfile ${pon} \
			--isnotusingponweights \
			--rawcov \
			--outputdir . \
			--log_directory .
	}	

	output {
		File seg = "cbs_seg/${sample_name}.seg"
		File copy_ratio = "tumor_pcov/${sample_name}.tn.tsv"
	}
}

task callSnvsOnPair {
	String tumor_name
	String normal_name
	File tumor_bam
	File normal_bam

	File tumor_bai = sub(tumor_bam, "\\.bam$", ".bai")
	File normal_bai = sub(normal_bam, "\\.bam$", ".bai")
	
	command {
		mutect \
			--analysis_type 'MuTect' \
			--intervals '/xchip/cga/reference/hg19/whole_exome_illumina_coding_v1_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list' \
			--tumor_sample_name ${tumor_name} \
			-I:tumor ${tumor_bam} \
			--normal_sample_name ${normal_name} \
			-I:normal ${normal_bam} \
			--reference_sequence '/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta' \
			--dbsnp '/xchip/cga/reference/hg19/dbsnp_134_b37.leftAligned.vcf' \
			--cosmic '/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf' \
			--normal_panel '/xchip/cga/reference/hg19/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf' \
			--out ${tumor_name}.call_stats.txt \
			--coverage_file ${tumor_name}.coverage.wig.txt \
			--power_file ${tumor_name}.power.wig.txt \
			--downsample_to_coverage '100000' \
			--enable_extended_output \
			--fraction_contamination '0.001'
	}

	output {
		File call_stats = "${tumor_name}.call_stats.txt"
		File coverage = "${tumor_name}.coverage.wig.txt"
		File power = "${tumor_name}.power.wig.txt"
	}
}

# select tumor SNVs that occur at normal heterozygous sites
task selectNormalHetSites {
	File in
	String out_fname
	File ref_snp

	command {
		callstats_ghet-filter.r \
			${in} \
			${out_fname} \
			--snp ${ref_snp}
	}

	output {
		File out = "${out_fname}"
	}
}

task convertCallStatsToCov {
	File in
	String out_fname
	
	command {
		covstats2cov.r \
			${in} \
			--output ${out_fname}
	}

	output {
		File out = "${out_fname}"
	}
}

task runAllelicCapseg {
	String sample_name
	File copy_ratio
	File seg
	File cov
	
	command {
		acapseg.r \
			--SID ${sample_name} \
			--capseg.probe.fn ${copy_ratio} \
			--capseg.seg.fn ${seg} \
			--germline.het.fn ${cov}
	}

	output {
		File out_rds = "${sample_name}.rds"
		File out_tsv = "results/${sample_name}.tsv"
	}
}
