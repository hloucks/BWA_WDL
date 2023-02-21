version 1.0

workflow MapReads{
    input {
        File fastq
        File reference
        File reference_bwt
        File reference_amb
        File reference_ann
        File reference_pac
        File reference_sa

        String ref_name=basename(sub(reference_amb, "\\.amb", ""))
        String fastq_name=basename(sub(sub(sub(fastq, "\\.gz$", ""), "\\.fastq$", ""), "\\.fa$", "")) #remove the file extension
        
        Int threadCount=20
        Int preemptible = 1
    }

    call FastqToBam {
        input:
            fastq=fastq,
            reference=reference,
            
            reference_bwt=reference_bwt,
            reference_amb=reference_amb,
            reference_ann=reference_ann,
            reference_pac=reference_pac,
            reference_sa=reference_sa,

            fastq_name=fastq_name,
            ref_name=ref_name,

            preemptible=preemptible,
            threadCount=threadCount
    }

    call SortBam {
        input:
            BamFile = FastqToBam.BamFile,
            fastq_name=fastq_name,
            ref_name=ref_name,

            preemptible=preemptible

    }


    output {
        File sortedBAM = SortBam.SortedBam
        File BamIndex = SortBam.BamIndex
    }

    parameter_meta {
        fastq: " Fastq reads to map to CHM13"
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }
}

task FastqToBam {
    input {
        File fastq
        File reference
        File reference_bwt
        File reference_amb
        File reference_ann
        File reference_pac
        File reference_sa

        String fastq_name
        String ref_name

        Int memSizeGB = 32
        Int preemptible
        Int threadCount
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # put my ref database files in one place 
        ln -s ~{reference}
        ln -s ~{reference_bwt}
        ln -s ~{reference_amb}
        ln -s ~{reference_ann}
        ln -s ~{reference_pac}
        ln -s ~{reference_sa}

        bwa mem -t ~{threadCount}  -c 128 -L 12,12 ~{ref_name} ~{fastq} | samtools view -Sb -@ 2 -O BAM -o ~{fastq_name}.~{ref_name}.bam - 

    >>>

    output {
        File BamFile="~{fastq_name}.~{ref_name}.bam"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/hdc-workflows/bwa-samtools"
    }
}

task SortBam {
    input {
        File BamFile

        String fastq_name
        String ref_name

        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools sort ~{BamFile} > ~{fastq_name}.~{ref_name}.sorted.bam
        samtools index ~{fastq_name}.~{ref_name}.sorted.bam

        
    >>>

    output {
        File SortedBam="~{fastq_name}.~{ref_name}.sorted.bam"
        File BamIndex="~{fastq_name}.~{ref_name}.sorted.bam.bai"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/hdc-workflows/bwa-samtools"
    }
}