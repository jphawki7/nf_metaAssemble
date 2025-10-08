#!/usr/bin/env nextflow

//Default Parameters
//For Read input - Can set data location to simplify input using id, or indicate full path to reads using analysis
params.data_location = "<YOUR_DATA_LOCATION>"
params.id = "reads"
params.analysis = "${params.data_location}/${params.id}"
params.reads = "${params.analysis}/*.fastq"
params.outdir = "${params.analysis}/analysis_output"
params.taxid = "1"
params.single = false
params.blastn = false
params.skip_diamond = false
params.skip_chop = false
params.skip_kraken = false
params.krakendb = "<YOUR_KRAKEN_DB>"
params.blastdb = "<YOUR_BLAST_DB>"
params.diamonddb = "<YOUR_DIAMOND_DB>"
params.hostile = false
params.help = false

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:

        nextflow run longboi.nf --id <subfolder in data_location with fastq reads> OR nextflow run longboi.nf --analysis <abs path to fastq reads>

       Typical arguments:
         --id                           Location of folder with fastq reads within data folder
         --outdir                       Location of output results (Default <analysis_folder>/analysis_output)
         --taxid                        Set taxid to extract reads for (Defaults 1 - Root - Will skip filtering)
       
       Toggle arguments:
         --skip_chop                    Skip Porechop and process reads through FASTP
         --single                       Set to only get requested taxid reads (Skips unclassified and root reads)
         --skip_diamond                 Disable DIAMOND search
         --skip_kraken                  Disable Kraken2 classification and read filtering
         --blastn                       Enable blastN search
         --hostile                      Abs path to host sequence for Hostile host removal - enables hostile
         
       Optional arguments:
         --analysis                     Abs path to folder with reads for analysis (Overides id) 
         --krakendb                     Abs path to krakenDB
         --blastdb                      Abs path to blastN database, used with --blastn
         --diamonddb                    Abs path to diamond database
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
    L O N G B O I - N F   P I P E L I N E
    ===================================
    Reads Dir    : ${params.reads}
    TaxID        : ${params.taxid}
    Kraken DB    : ${params.krakendb}
    outdir       : ${params.outdir}
    Diamond DB   : ${params.diamonddb}

    """
    .stripIndent(true)

process FASTQC {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads

    output:
    path "multiqc/*"

    script:
    """
    mkdir -p fastqc
    mkdir -p multiqc
    fastqc $reads -o fastqc -t ${task.cpus} --mem 2048 --nano
    multiqc fastqc/* -o multiqc
    """
}

process PORECHOP_ABI {
    tag "Processing ${reads.baseName}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads

    output:
    path "trimmed_reads/${reads.baseName}.fastq" //sets using basename command to generate .fq file

    script:
    """
    mkdir -p trimmed_reads
    porechop_abi -abi -i ${reads.baseName}.fastq -o trimmed_reads/${reads.baseName}.fastq -t ${task.cpus}
    """
}

process FASTP {
    tag "Quality trimming $reads"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads

    output:
    path "fastp/${reads.baseName}.fastq" //sets using basename command to generate .fq file

    script:
    """
    mkdir -p fastp
    fastplong -i ${reads.baseName}.fastq -o fastp/${reads.baseName}_original.fastq -w ${task.cpus} \
    -q 9 -l 50
    seqkit rmdup --ignore-case fastp/${reads.baseName}_original.fastq > fastp/${reads.baseName}.fastq
    """
}

process HOSTILE {
    tag "Dehosting $reads"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads

    output:
    path "hostile/${reads.baseName}.fastq" //sets using basename command to generate .fq file

    script:
    """
    mkdir -p hostile
    hostile clean --fastq1 ${reads.baseName}.fastq --out-dir hostile --index ${params.hostile}
    gunzip hostile/*.fastq.gz
    mv hostile/*.fastq hostile/${reads.baseName}.fastq
    """
}

process KRAKEN {
    tag "Running Kraken on $reads"
    publishDir params.outdir, mode: 'copy'
    label 'process_high_memory'
    
    input:
    path reads

    output:
    path "kraken2/*_out.tsv", emit: kraken_out 
    path "kraken2/*_report.tsv", emit: kraken_report
    path "kraken2/krona_${reads.baseName}.krona", emit: krona_input

    script:
    """
    mkdir -p kraken2
    kraken2 --threads ${task.cpus} --output kraken2/${reads.baseName}_out.tsv --report kraken2/${reads.baseName}_report.tsv --db ${params.krakendb} trimmed --confidence 0.1 $reads;
    kreport2krona.py -r kraken2/${reads.baseName}_report.tsv -o kraken2/krona_${reads.baseName}.krona
    """
}

process KRONA {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path report

    output:
    path "kraken2/krona_report.html"

    script:
    """
    mkdir -p kraken2
    ktImportText $report -o kraken2/krona_report.html
    """
}

process KRAKENFILTER {
    tag "Filtering $reads for ${params.taxid}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads
    path krakenout
    path kraken_report

    output:
    path "filtered_reads/*.fastq"
    path "combined_reads/*.fastq", emit: combined_reads

    script:
    def prefix = reads.baseName
    
    // Find the matching Kraken output and report files
    def matchedKrakenOut = krakenout.find { it.name.contains(prefix) }
    def matchedKrakenReport = kraken_report.find { it.name.contains(prefix) }
    
    if (!matchedKrakenOut || !matchedKrakenReport) {
        throw new RuntimeException("Matching Kraken2 files not found for prefix: ${prefix}")
    }
    """
    mkdir -p filtered_reads
    extract_kraken_reads.py -k ${matchedKrakenOut} -s $reads -o filtered_reads/${reads.baseName}_taxid_${params.taxid}.fastq --taxid ${params.taxid} --include-children -r ${matchedKrakenReport} --fastq-output;
    
    if [ ${params.single} = true ]; then
        echo "skipping root and unclassified reads"
    else
        echo "Extracting Root and Unclassified Reads"
        extract_kraken_reads.py -k ${matchedKrakenOut} -s $reads -o filtered_reads/${reads.baseName}_unclassified.fastq --taxid 0 -r ${matchedKrakenReport} --fastq-output;
        extract_kraken_reads.py -k ${matchedKrakenOut} -s $reads -o filtered_reads/${reads.baseName}_root.fastq --taxid 1 -r ${matchedKrakenReport} --fastq-output;
    fi

    mkdir -p combined_reads
    cat filtered_reads/${prefix}_*.fastq > combined_reads/${prefix}.fastq
    """
}

process METAMDBG {
    tag "metaMDBG assembly of ${reads.baseName}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads

    output:
    path "assemblies/${reads.baseName}/${reads.baseName}_metamdbg.fasta"

    script:
    """
    mkdir -p assemblies/${reads.baseName}
    metaMDBG asm --out-dir assemblies/${reads.baseName} --in-ont $reads --threads ${task.cpus};
    gunzip assemblies/${reads.baseName}/contigs.fasta.gz
    mv assemblies/${reads.baseName}/contigs.fasta assemblies/${reads.baseName}/${reads.baseName}_metamdbg.fasta
    """
}

process CANU {
    tag "canu assembly of ${reads.baseName}"
    publishDir params.outdir, mode: 'copy'
    label 'error_ignore'

    input:
    path reads

    output:
    path "assemblies/${reads.baseName}/${reads.baseName}_canu.fasta", emit: canu_contigs
    path "assemblies/${reads.baseName}/${reads.baseName}.unassembled_canu.fasta"

    script:
    """
    mkdir -p assemblies/${reads.baseName}
    canu -d assemblies/${reads.baseName} -p ${reads.baseName} -useGrid=false corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200 correctedErrorRate=0.105 genomeSize=0.01m corMaxEvidenceCoverageLocal=10 corMaxEvidenceCoverageGlobal=10 -trimmed -nanopore $reads
    mv assemblies/${reads.baseName}/${reads.baseName}.contigs.fasta assemblies/${reads.baseName}/${reads.baseName}_canu.fasta
    mv assemblies/${reads.baseName}/${reads.baseName}.unassembled.fasta assemblies/${reads.baseName}/${reads.baseName}.unassembled_canu.fasta
    #Remove Canu File if assembly is empty (Prevents Blasting Errors)
    [ -s "assemblies/${reads.baseName}/${reads.baseName}_canu.fasta" ] || rm "assemblies/${reads.baseName}/${reads.baseName}_canu.fasta"
    """
}

process FLYE {
    tag "Flye assembly of ${reads.baseName}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path reads

    output:
    path "assemblies/${reads.baseName}/${reads.baseName}_flye.fasta"

    script:
    """
    mkdir -p assemblies/${reads.baseName}
    flye --nano-hq $reads -o assemblies/${reads.baseName} --meta -t ${task.cpus}
    mv assemblies/${reads.baseName}/assembly.fasta assemblies/${reads.baseName}/${reads.baseName}_flye.fasta || true
    """
}

process DIAMOND_METAMDBG {
    tag "Blastx of ${contigs.baseName}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path contigs

    output:
    path "diamond/${contigs.baseName.split('_')[0]}/*.tsv"

    script:
    def filename = contigs.baseName.split('_')[0]

    """
    mkdir -p diamond/$filename;
    diamond blastx -d ${params.diamonddb} -q $contigs -o diamond/$filename/${contigs.baseName}.tsv \\
     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle --evalue 1e-6 -k 3 \\
     --threads ${task.cpus};
    #Insert header into blast file#
    sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle' diamond/$filename/${contigs.baseName}.tsv
    """
}

process DIAMOND_CANU {
    tag "Blastx of ${contigs.baseName}"
    publishDir params.outdir, mode: 'copy'
    label 'error_ignore'
    
    input:
    path contigs

    output:
    path "diamond/${contigs.baseName.split('_')[0]}/*.tsv"

    script:
    def filename = contigs.baseName.split('_')[0]

    """
    mkdir -p diamond/$filename;
    diamond blastx -d ${params.diamonddb} -q $contigs -o diamond/$filename/${contigs.baseName}.tsv \\
     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle --evalue 1e-6 -k 3 \\
     --threads ${task.cpus};
    #Insert header into blast file#
    sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle' diamond/$filename/${contigs.baseName}.tsv
    """
}
process DIAMOND_FLYE {
    tag "Blastx of ${contigs.baseName}"
    publishDir params.outdir, mode: 'copy'
    label 'error_ignore'
    
    input:
    path contigs

    output:
    path "diamond/${contigs.baseName.split('_')[0]}/*.tsv"

    script:
    def filename = contigs.baseName.split('_')[0]

    """
    mkdir -p diamond/$filename;
    diamond blastx -d ${params.diamonddb} -q $contigs -o diamond/$filename/${contigs.baseName}.tsv \\
     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle --evalue 1e-6 -k 3 \\
     --threads ${task.cpus};
    #Insert header into blast file#
    sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle' diamond/$filename/${contigs.baseName}.tsv
    """
}

process BLASTN_METAMDBG {
    tag "blasting ${contigs.baseName}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path contigs

    output:
    path "blast/${contigs.baseName.split('_')[0]}/*.tsv"

    script:
    def filename = contigs.baseName.split('_')[0]

    """
    mkdir -p blast/$filename;
    blastn -db ${params.blastdb} -query $contigs \\
     -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxid ssciname" \\
     -evalue 1e-6 -max_target_seqs 1 -out blast/$filename/${contigs.baseName}.tsv -num_threads ${task.cpus};
    #Insert header into blast file#
    sed -i '1i qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tstaxid\tssciname' blast/$filename/${contigs.baseName}.tsv
    """
}

process BLASTN_CANU {
    tag "blasting ${contigs.baseName}"
    publishDir params.outdir, mode: 'copy'
    label 'error_ignore'
    
    input:
    path contigs

    output:
    path "blast/${contigs.baseName.split('_')[0]}/*.tsv"

    script:
    def filename = contigs.baseName.split('_')[0]

    """
    mkdir -p blast/$filename;
    blastn -db ${params.blastdb} -query $contigs \\
     -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxid ssciname" \\
     -evalue 1e-6 -max_target_seqs 1 -out blast/$filename/${contigs.baseName}.tsv -num_threads ${task.cpus};    #Insert header into blast file#
    sed -i '1i qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tstaxid\tssciname' blast/$filename/${contigs.baseName}.tsv
    """
}

process BLASTN_FLYE{
    tag "blasting ${contigs.baseName}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path contigs

    output:
    path "blast/${contigs.baseName.split('_')[0]}/*.tsv"

    script:
    def filename = contigs.baseName.split('_')[0]

    """
    mkdir -p blast/$filename;
    blastn -db ${params.blastdb} -query $contigs \\
     -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxid ssciname" \\
     -evalue 1e-6 -max_target_seqs 1 -out blast/$filename/${contigs.baseName}.tsv -num_threads ${task.cpus};    #Insert header into blast file#
    sed -i '1i qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tstaxid\tssciname' blast/$filename/${contigs.baseName}.tsv
    """
}

workflow {

//QC Steps
    Channel
    .fromPath(params.reads)
    .collect()
    .set{inputFASTQC_ch}

    FASTQC(inputFASTQC_ch)
    
//Main Workflow
    Channel
    .fromPath(params.reads)
    .set{inputREADS_ch}

    if (!params.skip_chop) {
    porechop_ch = PORECHOP_ABI(inputREADS_ch)
    fastp_ch = FASTP(porechop_ch)
    } else {
    fastp_ch = FASTP(inputREADS_ch)
    }
    
    if (params.hostile){
        processed_reads_ch = HOSTILE(fastp_ch)
    } else {
        processed_reads_ch = fastp_ch
    }

    if (params.skip_kraken) {
        assembly_ch = processed_reads_ch
    } else{
        kraken_ch = KRAKEN(processed_reads_ch)
        krakenOut_ch = KRAKEN.out.kraken_out.collect()
        krakenReport_ch = KRAKEN.out.kraken_report.collect()
        krakenKrona_ch = KRAKEN.out.krona_input.collect()
        KRONA(krakenKrona_ch)
//Toggle for read filtering if taxid is changed from 1, only if kraken2 is enabled
            if (params.taxid != "1") {
                KRAKENFILTER(processed_reads_ch,krakenOut_ch,krakenReport_ch)
                assembly_ch = KRAKENFILTER.out.combined_reads
            } else {
                assembly_ch = processed_reads_ch
    }
    }

    metamdbg_ch = METAMDBG(assembly_ch)
    canu_ch = CANU(assembly_ch)
    canu_contigs_ch = CANU.out.canu_contigs
    flye_ch = FLYE(assembly_ch)
    
    if (!params.skip_diamond) {
            DIAMOND_METAMDBG(metamdbg_ch)
            DIAMOND_CANU(canu_contigs_ch)
            DIAMOND_FLYE(flye_ch)
    }
            
    if (params.blastn) {
            BLASTN_METAMDBG(metamdbg_ch)
            BLASTN_CANU(canu_contigs_ch)
            BLASTN_FLYE(flye_ch)
    }
}

workflow.onComplete {

    log.info ( workflow.success ? "\nDone! Results can be found in --> $params.outdir\n" : "Oops .. something went wrong" )

}