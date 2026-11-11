version 1.0
workflow iDrop_ISSAAC_ATAC_V1 {
    input{
    File refATAC
    Array[Array[File]] Data_ATAC
    String Outdir = "."
    String Species = "mm10"
    String chrMT = "chrM"
    String genomeSize = "mm"
    File? Whitelist
    String? Readformat
    String? readStructure
    String? down_rate
    Boolean down = true
    Boolean tenx = true
    Boolean format = false
    Boolean mixspecies = false
    }
    
########################################


    String fai=refATAC+"/genome.fasta"
    String ref_index=refATAC+"genome.index"
    String tss=refATAC+"tss.bed"
    String promo=refATAC+"promoter.bed"
    String chrmt=chrMT
    String gtf=refATAC+"genes.gtf"
    String species=Species
    String chrome_size=refATAC+"chrom.sizes"
    String? defaultConfig = "newT1"
    String config = select_first([readStructure,defaultConfig])
    String dsci_whitelist = "/C4_scATAC_analysis/whitelist_V1_30bp.txt"
    String whitelist = select_first([Whitelist,dsci_whitelist])
    # String read_format = if config == "oldT1" then "bc:50:-1,r2:0:49" else if config == "500" then "bc:50:59,bc:66:-1,r2:0:49" else "bc:0:19,r1:20:-1"
    String dsci_read_format = "bc:0:29,r1:30:-1,r2:0:49"
    String read_format = select_first([Readformat,dsci_read_format])

########################################### 

    scatter(fq in Data_ATAC) {
    File atacfq1 = if tenx == true then fq[1] else fq[0]
    File atacfq2 = if tenx == true then fq[2] else fq[1]
    File atacbc = fq[0]
}

    Array[String] Fastq1_atac = select_all(atacfq1)
    Array[String] Fastq2_atac = select_all(atacfq2)
    Array[String] barcode_atac = select_all(atacbc)
    Array[String] barcode = if config == "oldT1" then Fastq2_atac else if config == "500" then Fastq2_atac else barcode_atac

    call chromap {
    input:
        fastq1=Fastq1_atac,
        fastq2=Fastq2_atac,
        outdir=Outdir,
        ref_index=ref_index,
        ref=fai,
        whitelist=whitelist,
        barcode=barcode,
        tenx=tenx,
        down=down,
        format=format,
        down_rate=down_rate,
        read_format=read_format
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
}
    ###### script  ######
task chromap{
  input{
  Array[String] fastq1
  Array[String] fastq2
  String outdir
  String ref_index
  String ref
  String whitelist
  String ?down_rate
  Boolean ?tenx
  Boolean ?format
  Boolean ?down
  String read_format
  Array[String] barcode
  Int cpp=20
  Int mem=100
  }
  
command <<<

   mkdir -p outs
   faSize -detailed ~{ref} | sort -k1,1 -k2,2n > outs/genome.sizes
if (~{down})
then
  seqtk sample -s 42 ~{sep=',' fastq1} ~{down_rate} | gzip > outs/sub_read_1.fq.gz
  seqtk sample -s 42 ~{sep=',' fastq2} ~{down_rate} | gzip > outs/sub_read_2.fq.gz

  if (~{tenx})
  then
    seqtk sample -s 42 ~{sep=',' barcode} ~{down_rate} | gzip > outs/sub_I2.fq.gz
  fi
fi

if (~{format})
then
   if (~{down})
   then
     chromap -t ~{cpp} --preset atac -x ~{ref_index} -r ~{ref} -1 outs/sub_read_1.fq.gz -2 outs/sub_read_2.fq.gz -b ~{true="outs/sub_I2.fq.gz" false="outs/sub_read_1.fq.gz" tenx } --barcode-whitelist ~{whitelist} --read-format ~{read_format} -o outs/fragments.tsv 2> outs/alignment_report.tsv
   else
     chromap -t ~{cpp} --preset atac -x ~{ref_index} -r ~{ref} -1 ~{sep=',' fastq1} -2 ~{sep=',' fastq2} -b ~{sep=',' barcode} --barcode-whitelist ~{whitelist} --read-format ~{read_format} -o outs/fragments.tsv 2> outs/alignment_report.tsv
   fi
else
   if (~{down})
   then
     chromap -t ~{cpp} --preset atac -x ~{ref_index} -r ~{ref} -1 outs/sub_read_1.fq.gz -2 outs/sub_read_2.fq.gz -b ~{true="outs/sub_I2.fq.gz" false="outs/sub_read_1.fq.gz" tenx } --barcode-whitelist ~{whitelist} -o outs/fragments.tsv 2> outs/alignment_report.tsv
   else
     chromap -t ~{cpp} --preset atac -x ~{ref_index} -r ~{ref} -1 ~{sep=',' fastq1} -2 ~{sep=',' fastq2} -b ~{sep=',' barcode} --barcode-whitelist ~{whitelist} -o outs/fragments.tsv 2> outs/alignment_report.tsv
   fi
fi

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
  Int cpp=10
  Int mem=800
  }
  
command <<<

   mkdir -p outs
   mkdir -p macs2_pk
   mkdir -p outs/raw_mtx
   mkdir -p outs/filtered_mtx
   mkdir -p tmp
   mkdir -p temp
   cp ~{fragments} outs/fragments.tsv.gz
   python /ISSAAC/scripts/convert_frag_to_reads.py > tmp/temp.unsorted.bed
   bedClip tmp/temp.unsorted.bed ~{genomesize} tmp/temp.clipped.bed
   sort -k1,1 -k2,2n -T temp tmp/temp.clipped.bed | gzip > outs/reads.bed.gz
   rm -r tmp
   rm -r temp
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
  Int cpp=3
  Int mem=50
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
  Int mem=50
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
  Int mem=30
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
  Int mem=30
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
  Int cpp=3
  Int mem=100
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
  Int cpp=3
  Int mem=100
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