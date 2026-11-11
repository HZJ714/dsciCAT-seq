version 1.0
workflow iDrop_ISSAAC_full_10x_V1 {
    input{
    File refRNA
    File refATAC
    Array[Array[File]] Data_ATAC
    Array[Array[File]] Data_RNA
    String Outdir = "."
    String Species = "GRCh38-mm10"
    String chrMT = "GRCh38_chrM,mm10_chrM"
    String genomeSize = "5821961402"
    File Whitelist
    Boolean mixspecies = true
    }
    
########################################

    String Refdir=refRNA
    String Gtf=refRNA+"genes.gtf"
    String fai=refATAC+"/genome.fasta"
    String ref_index=refATAC+"genome.index"
    String tss=refATAC+"tss.bed"
    String promo=refATAC+"promoter.bed"
    String chrmt=chrMT
    String gtf=refATAC+"genes.gtf"
    String species=Species
    String chrome_size=refATAC+"chrom.sizes"
    String whitelist=Whitelist

########################################### 

    scatter(fq in Data_ATAC){
            File atacfq1 = fq[1]
            File atacfq2 = fq[2]
            File atacbc = fq[0]
        }
    scatter(fq in Data_RNA){
            File rnafq1 = fq[1]
            File rnafq2 = fq[2]
            File rnabc = fq[0]
        }
    
    Array[String] cDNA_Fastq1 = select_all(rnafq1)
    Array[String] cDNA_Fastq2 = select_all(rnafq2)
    Array[String] cDNA_barcode = select_all(rnabc)
    Array[String] Fastq1_atac = select_all(atacfq1)
    Array[String] Fastq2_atac = select_all(atacfq2)
    Array[String] barcode_atac = select_all(atacbc)

        call chromap {
            input:
                fastq1 = Fastq1_atac,
                fastq2 = Fastq2_atac,
                barcode = barcode_atac,
                outdir = Outdir,
                ref_index = ref_index,
                ref = fai,
                whitelist = whitelist
        }

    call atac_common {
    input:
        outdir=Outdir,
        fragments=chromap.fragments,
        genomesize=chromap.genomesize,
        genomeSize=genomeSize,
        chrMT=chrMT
    }

    if (mixspecies==false) {
        call atac_singlespecies1 {
            input:
            fragments=chromap.fragments
        }
    }

    if (mixspecies==false) {
        call atac_singlespecies2 {
            input:
            mtx=atac_common.mtx,
            bc=atac_common.bc
        }
    }

        if (mixspecies==true) {
        call atac_mixspecies1 {
            input:
            fragments=chromap.fragments
        }
    }

        if (mixspecies==true) {
        call atac_mixspecies2 {
            input:
            fragments=chromap.fragments
        }
    }

        if (mixspecies==true) {
        call atac_mixspecies3 {
            input:
            fragments=chromap.fragments
        }
    }

        if (mixspecies==true) {
        call atac_mixspecies4 {
            input:
            fragments=chromap.fragments
        }
    }
    
    call prepare_cb_umi {
            input:
            fastq2=cDNA_Fastq2,
            bc=cDNA_barcode
    }

        call rna_common {
            input:
                fastq1 = cDNA_Fastq1,
                cb_umi = prepare_cb_umi.cb_umi,
                whitelist = whitelist,
                refdir = Refdir
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

    ###### script  #####
task chromap{
  input{
  Array[String] fastq1
  Array[String] fastq2
  String outdir
  String ref_index
  String ref
  String whitelist
  Array[String] barcode
  Int cpp=20
  Int mem=100
  }
  
command <<<

   mkdir -p outs
   faSize -detailed ~{ref} | sort -k1,1 -k2,2n > outs/genome.sizes
   chromap -t ~{cpp} --preset atac -x ~{ref_index} -r ~{ref} -1 ~{sep=',' fastq1} -2 ~{sep=',' fastq2} -b ~{sep=',' barcode} --barcode-whitelist ~{whitelist} -o outs/fragments.tsv 2> outs/alignment_report.tsv
   bgzip outs/fragments.tsv && tabix -p bed outs/fragments.tsv.gz
   
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File fragments = "outs/fragments.tsv.gz"
    File fragments_index = "outs/fragments.tsv.gz.tbi"
    File aln_report = "outs/alignment_report.tsv"
    File genomesize  = "outs/genome.sizes"
  }
}

task atac_common{
  input{
  String outdir
  String fragments
  String genomesize
  String genomeSize
  String chrMT
  Int cpp=20
  Int mem=300
  }
  
command <<<

   mkdir -p outs
   mkdir -p macs2_pk
   mkdir -p tmp
   mkdir -p outs/raw_mtx
   mkdir -p outs/filtered_mtx
   cp ~{fragments} outs/fragments.tsv.gz
   /ISSAAC/scripts/convert_frag_to_reads.py | bedClip stdin ~{genomesize} stdout | sort -k1,1 -k2,2n | gzip > outs/reads.bed.gz
   macs2 callpeak -t outs/reads.bed.gz -g ~{genomeSize} -f BED --nomodel --extsize 200 --keep-dup all -B --SPMR --outdir macs2_pk -n aggregate 2> macs2_pk/peakcall_report.tsv
   cut -f 1-4 macs2_pk/aggregate_peaks.narrowPeak | sed '/~{chrMT}/d' | sort -k1,1 -k2,2n > macs2_pk/aggregate_peaks_formatted.bed
   bedtools intersect -a macs2_pk/aggregate_peaks_formatted.bed -b outs/reads.bed.gz -wo -sorted -g ~{genomesize} | sort -k8,8 | bedtools groupby -g 8 -c 4 -o freqdesc | gzip > outs/peak_read_ov.tsv.gz
   cut -f 1-3 macs2_pk/aggregate_peaks_formatted.bed > outs/raw_mtx/peaks.bed
   zcat outs/peak_read_ov.tsv.gz | cut -f 1 > outs/raw_mtx/barcodes.tsv
   /ISSAAC/scripts/generate_csc_mtx.py
   awk 'BEGIN{{OFS="\t"}}{{print $1 "-" $2 "-" $3, $1 "-" $2 "-" $3}}' outs/raw_mtx/peaks.bed > outs/raw_mtx/features.tsv
   STAR --runMode soloCellFiltering outs/raw_mtx outs/filtered_mtx/ --soloCellFilter EmptyDrops_CR
   
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File pk_report = "macs2_pk/peakcall_report.tsv"
    File narrowpeak = "macs2_pk/aggregate_peaks.narrowPeak"
    File formatpeak = "macs2_pk/aggregate_peaks_formatted.bed"
    File ov = "outs/peak_read_ov.tsv.gz"
    File raw_mtx = "outs/raw_mtx"
    File filtered_mtx = "outs/filtered_mtx"
    File mtx = "outs/filtered_mtx/matrix.mtx"
    File bc = "outs/filtered_mtx/barcodes.tsv"
  }
}

task atac_singlespecies1{
  input{
  String fragments
  Int cpp=2
  Int mem=35
  }
  
command <<<
   
   mkdir -p outs
   echo -e "size\tcount" > outs/fragment_size_distribution.tsv
   zcat ~{fragments} | sed '/^#/d' | awk '{print $3-$2}' | sort | uniq -c | sort -b -k2,2n | awk 'BEGIN{OFS="\t"}{print $2, $1}' >> outs/fragment_size_distribution.tsv
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File distribution = "outs/fragment_size_distribution.tsv"
  }
}

task atac_singlespecies2{
  input{
  String mtx
  String bc
  Int cpp=2
  Int mem=35
  }
  
command <<<
   
   mkdir -p outs
   /ISSAAC/scripts/output_qc_mtx  ~{mtx} ~{bc} outs/metrics.csv
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File metrics = "outs/metrics.csv"
  }
}

task atac_mixspecies1{
  input{
  String fragments
  Int cpp=2
  Int mem=35
  }
  
command <<<

   mkdir -p outs
   echo -e "size\tcount" > outs/fragment_size_distribution_GRCh38.tsv
   zcat ~{fragments} | sed '/^#/d' | grep GRCh38 | awk '{{print $3-$2}}' | sort | uniq -c | sort -b -k2,2n | awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' >> outs/fragment_size_distribution_GRCh38.tsv

   
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File distribution_hg = "outs/fragment_size_distribution_GRCh38.tsv"
  }
}

task atac_mixspecies2{
  input{
  String fragments
  Int cpp=2
  Int mem=35
  }
  
command <<<

   mkdir -p outs
   echo -e "size\tcount" > outs/fragment_size_distribution_mm10.tsv
   zcat ~{fragments} | sed '/^#/d' | grep mm10 | awk '{{print $3-$2}}' | sort | uniq -c | sort -b -k2,2n | awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' >> outs/fragment_size_distribution_mm10.tsv
   
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File distribution_mm = "outs/fragment_size_distribution_mm10.tsv"
  }
}

task atac_mixspecies3{
  input{
  String fragments
  Int cpp=2
  Int mem=50
  }
  
command <<<

   mkdir -p outs
   zcat ~{fragments} | grep GRCh38 | cut -f 4 | sort | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}' > outs/per_cell_barcode_total_fragment_count_GRCh38.tsv
   
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File frag_hg = "outs/per_cell_barcode_total_fragment_count_GRCh38.tsv"
  }
}

task atac_mixspecies4{
  input{
  String fragments
  Int cpp=2
  Int mem=50
  }
  
command <<<

   mkdir -p outs
   zcat ~{fragments} | grep mm10 | cut -f 4 | sort | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}' > outs/per_cell_barcode_total_fragment_count_mm10.tsv
   
   >>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File frag_mm = "outs/per_cell_barcode_total_fragment_count_mm10.tsv"
  }
}

task prepare_cb_umi{
        input{
        Array[String] fastq2
        Array[String] bc
        Int cpu=5
        Int mem=50
        }
        command <<<
        mkdir -p outs
        zcat ~{sep=',' bc} > /tmp/bc.tmp
        zcat ~{sep=',' fastq2} > /tmp/fastq2.tmp
        paste /tmp/bc.tmp /tmp/fastq2.tmp | awk -F '\t' '{{if(NR%4==1||NR%4==3) {{print $1}} if(NR%4==2||NR%4==0) {{print $1 substr($2,1,10)}} }}' | gzip > outs/CB_UMI.fastq.gz
        rm /tmp/bc.tmp /tmp/fastq2.tmp

>>>
runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }

        output{
                File cb_umi = "outs/CB_UMI.fastq.gz"
        }
}

task rna_common{
        input{
        Array[String] fastq1
        String cb_umi
        String whitelist
        String refdir
        Int cpu=20
        Int mem=200
        }
        command <<<
        mkdir -p outs
        STAR --genomeDir ~{refdir} --outFileNamePrefix outs/ --limitBAMsortRAM 40000000000 --readFilesCommand zcat --readFilesIn ~{sep=',' fastq1} ~{cb_umi} --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloType CB_UMI_Simple --clip3pNbases 115 --soloFeatures Gene --soloCBwhitelist ~{whitelist} --soloCellFilter EmptyDrops_CR --runThreadN ~{cpu} --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate

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
        Int cpu=2
        Int mem=35
        }
        command <<<
        mkdir -p outs
        samtools view -q 30 -F 3844 ~{bam} | grep GRCh38 | egrep -o 'CB:Z:[ACGT]{16}' | cut -f 3 -d ':' | sort | uniq -c  | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}' > outs/per_cell_barcode_total_read_count_GRCh38.tsv

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
        Int cpu=2
        Int mem=35
        }
        command <<<
        mkdir -p outs
        samtools view -q 30 -F 3844 ~{bam} | grep mm10 | egrep -o 'CB:Z:[ACGT]{16}' | cut -f 3 -d ':' | sort | uniq -c  | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}' > outs/per_cell_barcode_total_read_count_mm10.tsv

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
        Int cpu=2
        Int mem=35
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