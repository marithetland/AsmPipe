nextflow.enable.dsl=2

//READS-CHANNEL + CHECK IF ANY READS
reads_ch = Channel
        .fromFilePairs([params.reads_type1, params.reads_type2], flat: true, size: -1).ifEmpty {
        exit 1, "ERROR: did not find any read files with './*.fastq'"
        }


//CHECK IF BOTH READ1 AND READ2 ARE PROVIDED
reads_check_ch = reads_ch.map {
        if (it.size() != 3) {
                exit 1, "ERROR, didnt get exactly two readsets prefixed with ${it[0]}"
        }
        }

//Because of any(): True if any record, False if empty.
//Change to all(): False if any false record, True if empty.

process LOADFASTQ {

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        stdout

        script:
        """
        check_fastq.py --fastq $reads1
        check_fastq.py --fastq $reads2
        """
}
//add if unicycler048 and if unicycler050
//VERSION.TXT
process VERSIONS {
        
        publishDir path:("QC"), mode: 'copy'

        output:
        path("versions.txt")

        script:
        """
        echo "Program\tVersion" >> versions.txt
        unicycler --version >> versions.txt
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
        trim_galore --paired $reads1 $reads2
        """
}


//UNICYCLER

//Only one unicycler-process
//In script: if (params.unicycler048) {
//$params.path/to/unicycler048 -1 $reads1.......
//else if (params.unicycler050) {
//$params.path/to/unicycler048 -1 $reads1.......
//else {
//exit 1, "ERROR"}
//(Conditional scripts, under processes)
process ASSEMBLY {

        errorStrategy "${params.failure_action}"
        publishDir path:("fasta"), mode: 'copy', pattern: '*.fasta'
        publishDir path:("gfa_unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_assembly.gfa"}, pattern: '*.gfa'
        publishDir path:("logs/unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_unicycler.log"}, pattern: 'unicycler/unicycler.log'

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: fasta_files
        path("${sample_id}_assembly.gfa")
        path("unicycler/unicycler.log")


        script:
        """
        unicycler -1 $reads1 -2 $reads2 -o unicycler --verbosity 2 --keep 2 --depth_filter $params.depth_filter
        mv unicycler/assembly.fasta ${sample_id}_assembly.fasta
        mv unicycler/assembly.gfa ${sample_id}_assembly.gfa
        """
}


//UNICYCLER048
process UNICYCLER048 {

        errorStrategy "${params.failure_action}"
        conda "$params.unicycler048_env"

        publishDir path:("fasta"), mode: 'copy', pattern: '*.fasta'
        publishDir path:("gfa_unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_assembly.gfa"}, pattern: '*.gfa'
        publishDir path:("logs/unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_unicycler.log"}, pattern: 'unicycler/unicycler.log'

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: fasta_files
        path("${sample_id}_assembly.gfa")
        path("unicycler/unicycler.log")


        script:
        """
        unicycler -1 $reads1 -2 $reads2 -o unicycler --verbosity 2 --keep 2 --depth_filter $params.depth_filter --spades_path 
        mv unicycler/assembly.fasta ${sample_id}_assembly.fasta
        mv unicycler/assembly.gfa ${sample_id}_assembly.gfa
        """
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

        publishDir path:("QC"), mode: 'copy', pattern: 'multiqc_report.html'
        publishDir path:("QC/reports"), mode: 'copy',  saveAs: {filename -> "multiqc_fastqc.txt"} ,pattern: 'multiqc_data/multiqc_fastqc.txt'

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

        publishDir path:("QC/reports"), mode: 'copy', pattern: 'fast_count.tsv'

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

        publishDir path:("QC/reports"), mode: 'copy', saveAs: {filename -> "quast_transposed_report.tsv"}, pattern: 'quast_results/results*/transposed_report.tsv'
        publishDir path:("QC"), mode: 'copy', saveAs: {filename -> "quast_report.html"}, pattern: 'quast_results/results*/report.html'

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

        publishDir path:("QC/reports"), mode: 'copy'

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
process PANDAS {

        errorStrategy "${params.failure_action}"
        publishDir path:("QC"), mode: 'copy'

        input:
        tuple path(mlst), path(quast), path(multiqc), path(fast_count)

        output:
        path("final_report.tsv")

        //take asmbl_pandas.py in as params in config?
        script:
        """
        asmbl_pandas.py --mlst $mlst --quast $quast --multiqc $multiqc --fast_count $fast_count
        """
}

process RMLST {

        errorStrategy "${params.failure_action}" 
        publishDir path:("QC/rMLST"), mode: 'copy'

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
        publishDir path:("QC/kmerfinder"), mode: 'copy',  saveAs: {filename -> "${sample_id}_kmerfinder.csv"}

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
        publishDir path:("QC/Kleborate"), mode: 'copy'

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

        //MAKE VERSIONS.TXT
        VERSIONS()

        //FASTQ-INPUT CHECK
        LOADFASTQ(reads_ch)

        //RUN TRIM_GALORE
        if ( params.trim ) {
                
                TRIMMING(reads_ch)
                trimmed_ch = TRIMMING.out.trimmed_files
        }
        //SKIP TRIM_GALORE
        else {
                trimmed_ch = reads_ch
        }

        //RUN UNICYCLER048
        if ( params.unicycler048 ) {
                UNICYCLER048(trimmed_ch)
                assembly_ch = UNICYCLER048.out.fasta_files
        }
        //RUN NEWEST UNICYCLER
        else {
                ASSEMBLY(trimmed_ch)
                assembly_ch = ASSEMBLY.out.fasta_files
        }

        //MAKE LIST OF FASTQ FOR FAST_COUNT
        fastq_list_ch = reads_ch.map { it.drop(1) }.collect()
        //RUN FAST_COUNT
        fc_report = FASTCOUNT(fastq_list_ch)

        
        
        //MAKE LIST OF TRIMMED FASTQ FOR FASTQC
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

        //RUN PANDAS
        report_ch = mlst_report.concat( quast_report, multiqc_report, fc_report ).collect()
        PANDAS(report_ch)

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
