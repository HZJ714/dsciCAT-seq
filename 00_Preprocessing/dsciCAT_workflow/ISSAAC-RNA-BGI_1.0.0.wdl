version 1.0
workflow iDrop_ISSAAC_RNA_V1 {
    input{
    File refRNA
    Array[Array[File]] Data_RNA
    String Outdir = "."
    String Species = "GRCh38-mm10"
    String chrMT = "GRCh38_chrM,mm10_chrM"
    String genomeSize = "5821961402"
    File? Whitelist
    String? down_rate
    Boolean down = true
    Boolean mixspecies = true
    }
    
########################################

    String Refdir=refRNA
    String Gtf=refRNA+"genes.gtf"
    String chrmt=chrMT
    String species=Species
    String dsci_whitelist = "/C4_scATAC_analysis/whitelist_V1_30bp.txt"
    String whitelist = select_first([Whitelist,dsci_whitelist])

########################################### 

    scatter(fq in Data_RNA){
            File rnafq1 = fq[0]
            File rnafq2 = fq[1]
        }
    
    Array[String] cDNA_Fastq1 = select_all(rnafq1)
    Array[String] cDNA_Fastq2 = select_all(rnafq2)


    call rna_common {
            input:
            fastq1=cDNA_Fastq1,
            fastq2=cDNA_Fastq2,
            whitelist=whitelist,
            down=down,
            down_rate=down_rate,
            refdir=Refdir
    }

    if (mixspecies==false) {
        call rna_singlespecies {
            input:
            mtx_rna=rna_common.mtx_rna,
            bc_rna=rna_common.bc_rna
    }
   }
   
    if (mixspecies==true) {
        call rna_mixspecies1 {
            input:
            bam=rna_common.bam
    }
   }

    if (mixspecies==true) {
        call rna_mixspecies2 {
            input:
            bam=rna_common.bam
    }
   }
}
    ###### script  ######
task rna_common{
        input{
        Array[String] fastq1
        Array[String] fastq2
        String ?down_rate
        Boolean ?down
        String whitelist
        String refdir
        Int cpu=20
        Int mem=200
        }
        command <<<
        mkdir -p outs
        if (~{down})
        then
          seqtk sample -s 42 ~{sep=',' fastq1} ~{down_rate} | gzip > outs/sub_read_1.fq.gz
          seqtk sample -s 42 ~{sep=',' fastq2} ~{down_rate} | gzip > outs/sub_read_2.fq.gz
        fi
        
        if (~{down})
        then
        STAR --genomeDir ~{refdir} --outFileNamePrefix outs/ --limitBAMsortRAM 40000000000 --readFilesCommand zcat --readFilesIn outs/sub_read_2.fq.gz outs/sub_read_1.fq.gz --soloCBstart 1 --soloCBlen 30 --soloUMIstart 31 --soloUMIlen 12 --soloBarcodeReadLength 0 --clip3pNbases 0 --soloFeatures Gene --soloType CB_UMI_Simple --soloCBwhitelist ~{whitelist} --soloCellFilter EmptyDrops_CR --runThreadN ~{cpu} --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate
        else
          STAR --genomeDir ~{refdir} --outFileNamePrefix outs/ --limitBAMsortRAM 40000000000 --readFilesCommand zcat --readFilesIn ~{sep=',' fastq2} ~{sep=',' fastq1} --soloCBstart 1 --soloCBlen 30 --soloUMIstart 31 --soloUMIlen 12 --soloBarcodeReadLength 0 --clip3pNbases 0 --soloFeatures Gene --soloType CB_UMI_Simple --soloCBwhitelist ~{whitelist} --soloCellFilter EmptyDrops_CR --runThreadN ~{cpu} --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate
        fi

>>>
runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }

        output{
                File solo = "outs/Solo.out"
                File bam = "outs/Aligned.sortedByCoord.out.bam"
                File mtx_rna = "outs/Solo.out/Gene/filtered/matrix.mtx"
                File bc_rna = "outs/Solo.out/Gene/filtered/barcodes.tsv"       
        }
}

task rna_mixspecies1{
        input{
        String bam
        Int cpu=3
        Int mem=70
        }
        command <<<
        mkdir -p outs
        samtools view -q 30 -F 3844 ~{bam} | grep GRCh38 | egrep -o 'CB:Z:[ACGT]{30}' | cut -f 3 -d ':' | sort | uniq -c  | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}' > outs/per_cell_barcode_total_read_count_GRCh38.tsv

>>>
runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }

        output{
                File count_hg = "outs/per_cell_barcode_total_read_count_GRCh38.tsv"
        }
}

task rna_mixspecies2{
        input{
        String bam
        Int cpu=3
        Int mem=70
        }
        command <<<
        mkdir -p outs
        samtools view -q 30 -F 3844 ~{bam} | grep mm10 | egrep -o 'CB:Z:[ACGT]{30}' | cut -f 3 -d ':' | sort | uniq -c  | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}' > outs/per_cell_barcode_total_read_count_mm10.tsv

>>>
runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }

        output{
                File count_mm = "outs/per_cell_barcode_total_read_count_mm10.tsv"
        }
}

task rna_singlespecies{
        input{
        String mtx_rna
        String bc_rna
        Int cpu=3
        Int mem=70
        }
        command {
       mkdir -p outs/Solo.out
       /ISSAAC/scripts/output_qc_mtx  ~{mtx_rna} ~{bc_rna} outs/Solo.out/metrics.csv
}
runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }

        output{
                File metrics_rna = "outs/Solo.out/metrics.csv"
        }
}