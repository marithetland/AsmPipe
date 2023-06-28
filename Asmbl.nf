nextflow.enable.dsl=2

//READS_CH + CHECK IF ANY READS
reads_ch = Channel
        .fromFilePairs([params.reads_type1, params.reads_type2], flat: true, size: -1).ifEmpty {
        exit 1, "ERROR: did not find any read files with ${params.reads_type1} or ${params.reads_type2}."
        }

//CHECK IF BOTH READ1 AND READ2 ARE PROVIDED
paired_reads_check_ch = reads_ch.map {
        if (it.size() != 3) {
                exit 1, "ERROR, didnt get exactly two readsets prefixed with ${it[0]}."
        }
        }

//CHECK FILES, EMPTY OR NOT FASTQ-FORMAT
process LOADFASTQ {

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), stdout

        script:
        """
        check_fastq.py --fastq $reads1
        check_fastq.py --fastq $reads2
        """
}

//VERSION.TXT
process VERSIONS {
        publishDir path:("reports"), mode: 'copy'

        input:
        val(unicycler)
        output:
        path("versions.txt")

        script:
        """
        echo "Program\tVersion" >> versions.txt
        $unicycler --version >> versions.txt
        spades.py --version >> versions.txt
        trim_galore --version | grep version | tr -d " " | sed "s/^/trim_galore\t/g" >> versions.txt
        cutadapt --version | sed "s/^/cutadapt\t/g"  >> versions.txt
        fastqc --version >> versions.txt
        multiqc --version >> versions.txt
        mlst --version >> versions.txt
        quast.py --version >> versions.txt
        kleborate --version >> versions.txt
        """
}

//RENAME
process RENAME {
        publishDir path:("fastq"), mode: 'copy'
        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

        script:
        if ("$reads1" != "${sample_id}_1.fastq.gz"){
                """
                mv $reads1 ${sample_id}_1.fastq.gz
                mv $reads2 ${sample_id}_2.fastq.gz
                """
        }
        else {
                """
                echo $reads1
                echo $reads2
                """
        }

}

//TRIM_GALORE
process	TRIMMING {
        errorStrategy "${params.failure_action}"
        publishDir path:("trimmed_reads"), mode: 'copy', pattern: '*fq.gz'
        publishDir path:("logs/trim_galore"), mode: 'copy', pattern: '*_trimming_report.txt'

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), path("${sample_id}*_val_1.fq.gz"), path("${sample_id}*_val_2.fq.gz"), emit: trimmed_files
        path("${reads1}_trimming_report.txt")
        path("${reads2}_trimming_report.txt")

        script:
        """
        trim_galore --paired $reads1 $reads2 --cores 4
        """
}

//UNICYCLER
process ASSEMBLY {
        maxForks = "${params.maxForks}"
        errorStrategy "${params.failure_action}"
        publishDir path:("fasta"), mode: 'copy', pattern: '*.fasta'
        publishDir path:("gfa"), mode: 'copy', saveAs: {filename -> "${sample_id}_assembly.gfa"}, pattern: '*.gfa'
        publishDir path:("logs/unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_unicycler.log"}, pattern: 'unicycler/unicycler.log'

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: fasta_files
        path("${sample_id}_assembly.gfa")
        path("unicycler/unicycler.log")

        script:
        if (params.unicycler048) {
                """
                $params.unicycler048_path -1 $reads1 -2 $reads2 -o unicycler --verbosity 2 --keep 2 --depth_filter $params.depth_filter --pilon_path $params.unicycler048_path
                mv unicycler/assembly.fasta ${sample_id}_assembly.fasta
                mv unicycler/assembly.gfa ${sample_id}_assembly.gfa
                """
        }
        else {
                """
                $params.unicycler050_path -1 $reads1 -2 $reads2 -o unicycler --verbosity 2 --keep 2 --depth_filter $params.depth_filter
                mv unicycler/assembly.fasta ${sample_id}_assembly.fasta
                mv unicycler/assembly.gfa ${sample_id}_assembly.gfa
                """
        }
}

//FASTQC
process FASTQC {
        errorStrategy "${params.failure_action}"

        input:
        path(rawreads)

        output:
        path("*fastqc.zip")

        script:
        """
        fastqc $rawreads 
        """
}

//MULTIQC
process MULTIQC {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports"), mode: 'copy', pattern: 'multiqc_report.html'
        publishDir path:("reports/extra"), mode: 'copy',  saveAs: {filename -> "multiqc_fastqc.txt"} ,pattern: 'multiqc_data/multiqc_fastqc.txt'

        input:
        path(fastqc_files)
        
        output:
        path("multiqc_report.html")
        path("multiqc_data/multiqc_fastqc.txt"), emit: report

        script:
        """
        multiqc $fastqc_files
        """
}

//FAST_COUNT
process FASTCOUNT {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports/extra"), mode: 'copy', pattern: 'fast_count.tsv'

        input:
        path(rawreads)

        output: 
        path('fast_count.tsv')

        script:
        """
        fast_count >> fast_count.tsv
        fast_count $rawreads >> fast_count.tsv
        """
}

//QUAST
process QUAST {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports/extra"), mode: 'copy', saveAs: {filename -> "quast_transposed_report.tsv"}, pattern: 'quast_results/results*/transposed_report.tsv'
        publishDir path:("reports"), mode: 'copy', saveAs: {filename -> "quast_report.html"}, pattern: 'quast_results/results*/report.html'

        input:
        path(fasta)

        output:
        path("quast_results/results*/transposed_report.tsv"), emit: rep
        path("quast_results/results*/report.html")

        script:
        """
        quast.py $fasta
        """
}

//MLST
process MLST {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports/extra"), mode: 'copy'

        input:
        path(fasta)

        output:
        path("mlst.tsv")

        script:
        """
        mlst $fasta > mlst.tsv
        """
}

//FINAL_REPORT
process FINAL_REPORT {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports"), mode: 'copy'

        input:
        tuple path(mlst), path(quast), path(multiqc), path(fast_count)

        output:
        path("final_report.tsv")

        script:
        """
        asmbl_pandas.py --mlst $mlst --quast $quast --multiqc $multiqc --fast_count $fast_count
        """
}

process RMLST {
        errorStrategy "${params.failure_action}" 
        publishDir path:("reports/rMLST"), mode: 'copy'

        input:
        path(fasta)

        output:
        path("rMLST.tsv")

        script:
        """
        rmlst_script.py -f $fasta >> rMLST.tsv
        """
}

process KMERFINDER {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports/kmerfinder"), mode: 'copy',  saveAs: {filename -> "${sample_id}_kmerfinder.csv"}

        input:
        tuple val(sample_id), path(fasta)

        output:
        path("output/results.txt")

        script:
        """
        kmerfinder.py -i $fasta -db $params.ATG -tax $params.tax
        """
}

process KLEBORATE {
        errorStrategy "${params.failure_action}"
        publishDir path:("reports/Kleborate"), mode: 'copy'

        input:
        path(fasta)

        output:
        path("Kleborate_results.txt")

        script:
        """
        kleborate --all -a $fasta
        """
}

workflow {
        
        fastq_ch = RENAME(reads_ch)

        //MAKE VERSIONS.TXT
        if ( params.unicycler048 ) {
                VERSIONS(params.unicycler048_path)
        }
        else {
                VERSIONS(params.unicycler050_path)
        }

        //FASTQ-INPUT CHECK
        check_fastq_ch = LOADFASTQ(reads_ch). map {
                if (it =~ "False") {
                        exit 1, "ERROR, one of the fastq files prefixed with ${it[0]} is not in accordance with fastq format (either corrupt or empty)."
                }
        }

        //MAKE FASTQ_LIST FOR FAST_COUNT
        fastq_list_ch = fastq_ch.map { it.drop(1) }.collect()
        //RUN FAST_COUNT
        fc_report = FASTCOUNT(fastq_list_ch)

        //RUN TRIM_GALORE
        if ( params.trim ) {
                TRIMMING(fastq_ch)
                trimmed_ch = TRIMMING.out.trimmed_files
        }
        //SKIP TRIM_GALORE
        else {
                trimmed_ch = fastq_ch
        }

        //RUN UNICYCLER
        ASSEMBLY(trimmed_ch)
        assembly_ch = ASSEMBLY.out.fasta_files

        //MAKE TRIMMED_LIST FOR FASTQC
        trimmed_list_ch = trimmed_ch.map { it.drop(1) }.collect()
        //RUN FASTQC
        fastqc_ch = FASTQC(trimmed_list_ch)
        //RUN MULTIQC
        MULTIQC(fastqc_ch)
        multiqc_report = MULTIQC.out.report

        //RUN QUAST
        QUAST(assembly_ch.map { it.drop(1)}.collect())
        quast_report = QUAST.out.rep

        //RUN MLST
        mlst_report = MLST(assembly_ch.map { it.drop(1)}.collect())

        //MAKE FINAL REPORT
        report_ch = mlst_report.concat( quast_report, multiqc_report, fc_report ).collect()
        FINAL_REPORT(report_ch)

        //rMLST
        if ( params.rmlst ) {
                RMLST(assembly_ch.map { it.drop(1)}.collect())
        }

        //KMERFINDER
        if ( params.kmerfinder ) {
                KMERFINDER(assembly_ch)
        }
        
        //KLEBORATE
        if ( params.kleborate ) {
                KLEBORATE(assembly_ch.map { it.drop(1)}.collect())
        }
}