nextflow.enable.dsl=2



//add checkIfExists: true? .fromFilePairs(params.reads, checkIfExists: true)
reads_ch = Channel
        .fromFilePairs([params.reads_type1, params.reads_type2], flat: true)


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

// //UNICYCLER
// process ASSEMBLY {

//         errorStrategy "${params.failure_action}"

//         publishDir path:("fasta/unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_assembly.fasta"}, pattern: 'unicycler/*.fasta'
//         publishDir path:("gfa_unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_assembly.gfa"}, pattern: 'unicycler/*.gfa'
//         publishDir path:("logs/unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_unicycler.log"}, pattern: 'unicycler/unicycler.log'

//         input:
//         tuple val(sample_id), path(reads1), path(reads2)

//         output:
//         tuple val(sample_id), path("unicycler/assembly.fasta"), emit: fasta_files
//         path("unicycler/assembly.gfa")
//         path("unicycler/unicycler.log")


//         script:
//         """
//         unicycler -1 $reads1 -2 $reads2 -o unicycler --verbosity 2 --keep 2 --depth_filter $params.depth_filter
//         """
// }



//UNICYCLER
process ASSEMBLY {

        errorStrategy "${params.failure_action}"

        publishDir path:("fasta"), mode: 'copy', pattern: 'unicycler/*.fasta'
        publishDir path:("gfa_unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_assembly.gfa"}, pattern: 'unicycler/*.gfa'
        publishDir path:("logs/unicycler"), mode: 'copy', saveAs: {filename -> "${sample_id}_unicycler.log"}, pattern: 'unicycler/unicycler.log'

        input:
        tuple val(sample_id), path(reads1), path(reads2)

        output:
        tuple val(sample_id), path("unicycler/${sample_id}_assembly.fasta"), emit: fasta_files
        path("unicycler/${sample_id}_assembly.gfa")
        path("unicycler/unicycler.log")


        script:
        """
        unicycler -1 $reads1 -2 $reads2 -o unicycler --verbosity 2 --keep 2 --depth_filter $params.depth_filter
        mv unicycler/assembly.fasta unicycler/${sample_id}_assembly.fasta
        mv unicycler/assembly.gfa unicycler/${sample_id}_assembly.gfa
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

//POLYPOLISH
process POLYPOLISH {
        
        errorStrategy "${params.failure_action}"
        
        conda "$params.polypolish_env"

        publishDir path:("fasta/polypolish"), mode: 'copy', pattern: '*_polypolish.fasta'

        input:
        tuple val(sample_id), path(illumina1), path(illumina2), path(fasta)

        output:
        tuple val(sample_id), path("${sample_id}_polypolish.fasta")

        //add seqkit sort --by-length --reverse contigs.fasta ??
        //or
        //add seqkit sort --by-length --reverse contigs.fasta | seqkit replace --pattern '.+' --replacement 'Contig_{nr}' ??
        //include polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam ??
        script:
        """
        bwa index $fasta
        bwa mem -t 16 -a $fasta $illumina1 > ${sample_id}_1.sam
        bwa mem -t 16 -a $fasta $illumina2 > ${sample_id}_2.sam
        polypolish $fasta ${sample_id}_1.sam ${sample_id}_2.sam > ${sample_id}_polypolish.fasta
        """
}

//POLCA
process POLCA {

        errorStrategy "${params.failure_action}"

        publishDir path:("fasta/polca"), mode: 'copy', pattern: '*_polypolish_polca.fasta'

        conda "$params.polca_env"

        input:
        tuple val(sample_id), path(illumina1), path(illumina2), path(polypolished)

        output:
        path("${sample_id}_polypolish_polca.fasta")

        //if the project directories are going to be deleted, the new name made by using mv can instead be done through publishDir with "filename"
        script:
        """
        polca.sh -a $polypolished -r "$illumina1 $illumina2" -t 16 -m 1G
        mv ${sample_id}_polypolish.fasta.PolcaCorrected.fa ${sample_id}_polypolish_polca.fasta
        """
}

//FINAL_REPORT
process PANDAS {

publishDir path:("QC"), mode: 'copy'

        input:
        tuple path(mlst), path(quast), path(multiqc), path(fast_count)

        output:
        path("final_report.tsv")

        //take asmbl_pandas.py in as params in config?
        script:
        """
        python ~/Programs/Asmbl-development-nf/asmbl_pandas.py --mlst $mlst --quast $quast --multiqc $multiqc --fast_count $fast_count
        """
}


workflow {

        //RUN TRIM_GALORE (or not)
        if ( params.trim ) {
                
                TRIMMING(reads_ch)
                trimmed_ch = TRIMMING.out.trimmed_files
        }
        else {
                trimmed_ch = reads_ch
        }

        //RUN UNICYCLER
        ASSEMBLY(trimmed_ch)
        assembly_ch = ASSEMBLY.out.fasta_files

        //RUN POLYPOLISH
        assembly_trimmed_ch = trimmed_ch.combine(assembly_ch, by: 0)


        if ( params.polypolish ) {
                polypolished_ch = POLYPOLISH(assembly_trimmed_ch)
                
        }

        //RUN POLCA
        polypolished_trimmed_ch = trimmed_ch.combine(polypolished_ch, by: 0)

        if ( params.polca ) {
                polca_ch = POLCA(polypolished_trimmed_ch)
        }

        //MAKE LIST OF FASTQ FOR FAST_COUNT
        fastq_list_ch = reads_ch.map { it.drop(1) }.collect()
        //RUN FAST_COUNT
        fc_report = FASTCOUNT(fastq_list_ch)

        
        
        //MAKE LIST OF VAL-FASTQ FOR FASTQC
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

}
