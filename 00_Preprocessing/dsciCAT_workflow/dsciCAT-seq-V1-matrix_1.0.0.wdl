version 1.0
workflow iDrop_dsciCAT_1 {
    input{
    File refRNA
    File refATAC
    Array[Array[File]] Data_ATAC
    Array[Array[File]] Data_RNA
    String Outdir = "."
    String Species = "mm10"
    String chrMT = "chrM"
    String genomeSize = "mm"
    String ProjectID = "<project-id>"
    String SampleID
    File BeadsBarcode
    Int? ForceFrag
    Int? expectCellNum
    Int? forceCellNum
    String? readStructure
    Boolean Intron = true
    Boolean mixspecies = false
    }
    Int? minUMI = 1000
    String Root = "/dsciRNA"
    String ChrSpecies = "/dsciRNA/species_binding.txt"
    
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
    String macs2sp=genomeSize
    String root_atac = "/C4_scATAC_analysis"
    String ID = SampleID
    String PID = "sentieon.q -P "+ProjectID+"_sentieon"
    String outdir = Outdir
    String runID = ID
    String? defaultConfig = "newT1"
    String config = select_first([readStructure,defaultConfig])
    String whitelist = "/C4_scATAC_analysis/whitelist_V1_30bp.txt"
    # String read_format = if config == "oldT1" then "bc:50:-1,r2:0:49" else if config == "500" then "bc:50:59,bc:66:-1,r2:0:49" else "bc:0:19,r1:20:-1"
    String read_format = "bc:0:29,r1:30:-1,r2:0:49"

########################################### 

    scatter(fq in Data_ATAC){
            File atacfq1 = fq[0]
            File atacfq2 = fq[1]
        }
    scatter(fq in Data_RNA){
            File rnafq1 = fq[0]
            File rnafq2 = fq[1]
        }
    
    Array[String] cDNA_Fastq1 = select_all(rnafq1)
    Array[String] cDNA_Fastq2 = select_all(rnafq2)
    Array[String] Fastq1_atac = select_all(atacfq1)
    Array[String] Fastq2_atac = select_all(atacfq2)
    Array[String] barcode = if config == "oldT1" then Fastq2_atac else if config == "500" then Fastq2_atac else Fastq1_atac

    call mapping {
    input:
        fastq1=Fastq1_atac,
        fastq2=Fastq2_atac,
        outdir=Outdir,
        root_atac=root_atac,
        ref_index=ref_index,
        ref=fai,
        whitelist=whitelist,
        barcode=barcode,
        read_format=read_format
    }

    call deconvolution { 
    input: 
        alnbed=mapping.alnbed,
        outdir=Outdir,
        root_atac=root_atac,
        ID=ID,
        mixspecies = mixspecies,
        species = species,
        ForceFrag=ForceFrag,
        chrome_size=chrome_size,
        tss=tss,
        chrmt=chrmt 
    }
    call qc {
    input:
        root_atac=root_atac,
        tss=tss,
        FragmentFile=deconvolution.FragmentFile,
        outdir=outdir,
        ID=ID,
        chrmt=chrmt,
    }

    call peakcount {
    input:
        root_atac=root_atac,
        listtxt=deconvolution.listtxt,
        macs2sp= macs2sp,
        promo= promo,
        chrmt=chrmt,
        outdir=outdir,
        FragmentFile=deconvolution.FragmentFile,
        ID=ID,
        Metadata=deconvolution.Metadata
    }
    call cDNAanno{
            input:
            outdir=Outdir,
            root=Root,
            fastq1=cDNA_Fastq1,
            fastq2=cDNA_Fastq2,
            refdir=Refdir,
            chrSpecies=ChrSpecies,
            barcode=BeadsBarcode,
            barcode_counts_raw="cDNA_barcode_counts_raw.txt",
            report="cDNA_sequencing_report.csv",
            gtf=gtf,
            intron=Intron
    }
        call getM280UMI{
            input:
            outdir=Outdir,
            root=Root,
            SampleID=SampleID,
            finalsortbam=cDNAanno.finalsortbam,
            beads_stat=cDNAanno.beads_stat,
            rawmatrix=cDNAanno.rawmatrix,
            ATAC_barcodeTran=deconvolution.barcodeMerge
    }
    call count{
        input:
        outdir=Outdir,
        root=Root,
        anno_decon=getM280UMI.anno_decon,
        ID=ID,
        cellid=getM280UMI.cellid
    }
    if (Intron==true) {
        call splice_matrix{
            input:
            outdir=Outdir,
            root=Root,
            anno_decon=getM280UMI.anno_decon,
            countmatrix=count.matrixdir,
            cellid=getM280UMI.cellid
        }
    }
    call statSaturation{
            input:
            outdir=Outdir,
            root=Root,
            anno_decon=getM280UMI.anno_decon,
            ID=ID,
            barcodeTranslate=deconvolution.barcodeMerge,
            beads_stat=cDNAanno.beads_stat
    }
    }
    ###### script  ######
task mapping {
  input{
  Array[String] fastq1
  Array[String] fastq2
  String outdir
  String ref_index
  String ref
  String root_atac
  String whitelist
  String read_format
  Array[String] barcode
  Int cpp=10
  Int mem=100
  }
  
command {

   export PATH="/usr/local/bin:$PATH"
   export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib64/:$LD_LIBRARY_PATH"
   mkdir -p ATAC
   mkdir -p ATAC/02.chromap
   chromap --preset atac --bc-error-threshold 0 --trim-adapters -x ${ref_index} -r ${ref} -1 ${sep=',' fastq1} -2 ${sep=',' fastq2} -o ATAC/02.chromap/aln.bed --barcode ${sep=',' barcode} --barcode-whitelist ${whitelist} --read-format ${read_format} -t ${cpp} 2> ATAC/02.chromap/alignment_report.tsv
   }
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File alnbed = "ATAC/02.chromap/aln.bed"
    File aln_report = "ATAC/02.chromap/alignment_report.tsv"
  }
}

task deconvolution {
  input{
  String outdir
  String alnbed
  String root_atac
  String species
  String chrmt
  String chrome_size
  String tss
  String ID
  Boolean mixspecies
  Int ?ForceFrag
  Int ?mapq
  Int cpp=10
  Int mem=150
  }
  
  command <<<

   export PATH="/usr/local/bin:$PATH"
   export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib64/:$LD_LIBRARY_PATH"
   mkdir -p ATAC/03.d2cfile
   if (~{mixspecies})
   then
   ~{root_atac}/d2c_v1.5.3/bin/d2c merge -i ~{alnbed} --mapq 30 --bf ~{default=0 ForceFrag} -r ~{species} --mix-species -o ATAC/03.d2cfile -c 10 --tn5 10  -n ~{ID} --mc ~{chrmt}  --sat --bt1 CB --bg ~{chrome_size} --ts ~{tss} --rank inflection
   Rscript '--args -I ATAC/03.d2cfile/~{ID}.Metadata.tsv -O ATAC' <hpc-mount>/<project-id>/<user>/shell/root_atac/script/sta_collision_rate.R
   perl ~{root_atac}/script/fishInWinter.fq.gz.pl  -bf table -ff table -bc 1 -fc 4 ATAC/03.d2cfile/Mouse_barcode_list.txt ATAC/03.d2cfile/~{ID}.fragments.tsv.gz > ATAC/03.d2cfile/~{ID}.Mouse.fragments.tsv
    perl ~{root_atac}/script/fishInWinter.fq.gz.pl  -bf table -ff table -bc 1 -fc 4 ATAC/03.d2cfile/Human_barcode_list.txt ATAC/03.d2cfile/~{ID}.fragments.tsv.gz > ATAC/03.d2cfile/~{ID}.Human.fragments.tsv
    bgzip ATAC/03.d2cfile/~{ID}.Mouse.fragments.tsv
    bgzip ATAC/03.d2cfile/~{ID}.Human.fragments.tsv
    tabix -p bed ATAC/03.d2cfile/~{ID}.Mouse.fragments.tsv.gz 
    tabix -p bed ATAC/03.d2cfile/~{ID}.Human.fragments.tsv.gz
    tail -n +2 ATAC/03.d2cfile/~{ID}.Metadata.tsv |awk '{print $1}' > ATAC/03.d2cfile/list.~{ID}.txt
   else
   ~{root_atac}/d2c_v1.5.3/bin/d2c merge -i ~{alnbed} --mapq 30 --bf ~{default=0 ForceFrag} -r ~{species} -o ATAC/03.d2cfile -c 10 --tn5 10  -n ~{ID} --mc ~{chrmt}  --sat --bt1 CB --bg ~{chrome_size} --ts ~{tss} --rank inflection
    cp ATAC/03.d2cfile/~{ID}.CorrelationBarcodes.tsv.gz ATAC/report/plot_input2_Jaccard_Overlap_Knee.csv.gz
    cp ATAC/03.d2cfile/~{ID}.barcodeCount.tsv ATAC/report/plot_input1_Bead_Barcode_Knee.csv
    tail -n +2 ATAC/03.d2cfile/~{ID}.Metadata.tsv |awk '{print $1}' > ATAC/03.d2cfile/list.~{ID}.txt
   fi

>>>
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File FragmentFile="ATAC/03.d2cfile/~{ID}.fragments.tsv.gz"
    File FragmentFileIndex="ATAC/03.d2cfile/~{ID}.fragments.tsv.gz.tbi"
    File listtxt="ATAC/03.d2cfile/list.~{ID}.txt"
    File Correlation="ATAC/03.d2cfile/~{ID}.CorrelationBarcodes.tsv.gz"
    File barcodeCount="ATAC/03.d2cfile/~{ID}.barcodeCount.tsv"
    File barcodeMerge="ATAC/03.d2cfile/~{ID}.barcodeMerge.tsv"
    File Metadata="ATAC/03.d2cfile/~{ID}.Metadata.tsv"
    File d2cCutoff="ATAC/03.d2cfile/~{ID}.d2cCutoff.tsv"
    File sequenceSaturation="ATAC/03.d2cfile/~{ID}.sequenceSaturation.tsv"
  }
}

task qc {
    input{
    String outdir
    String root_atac
    String tss
    String FragmentFile
    String ID
    String chrmt
    Int cpp=3
    Int mem=60
    }

  command <<<
        mkdir -p ATAC/report
        Rscript ~{root_atac}/script/plot_TSSEnrichment_FragSize_nochrM.R -T ~{tss} -F ~{FragmentFile} -G ~{ID} -O ATAC -MT ~{chrmt}
        nf=$(gunzip -c ~{FragmentFile} | awk '{if($3-$2<147) print $0}' | wc -l)
        mn=$(gunzip -c ~{FragmentFile} | awk '{if($3-$2>147 && $3-$2<294) print $0}' | wc -l)
        total=$(gunzip -c ~{FragmentFile} | wc -l)

        printf "$nf $total" | awk '{print "Fraction of nucleosome-free regions:"$1/$2*100"%"}' >> ATAC/report/5_2.library.QC.csv
        printf "$mn $total" | awk '{print "Fraction of fragments mono-nucleosome regions:"$1/$2*100"%"}' >> ATAC/report/5_2.library.QC.csv

  >>>
  runtime {
    docker_url:"<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File Qcfie="ATAC/report/5_2.library.QC.csv"
    File plot6_TSS="ATAC/report/plot6_TSS.svg"
    File plot5_InterSize="ATAC/report/plot5_InterSize.svg"
    
  }
}


task peakcount {
        input{
        String root_atac
        String outdir
        String macs2sp
        String promo
        String ID
        String FragmentFile
        String listtxt
        String chrmt
        String Metadata
        Int cpp=8
        Int mem=100
        }


  command <<<
    mkdir -p ATAC/out
    mkdir -p ATAC/out/Peak
    mkdir -p ATAC/out/Promoter
    macs2 callpeak -t ~{FragmentFile} -f BED -g ~{macs2sp} -n ~{ID} -B -q 0.001 --nomodel --outdir ATAC/out
    Rscript ~{root_atac}/script/C4scATAC_Cluster_AnnotationV3.R -I ATAC/out/~{ID}_peaks.narrowPeak -F ~{FragmentFile} -C ~{listtxt} -G ~{promo} -MT ~{chrmt} -Q ~{Metadata} -O ATAC
  >>>
  runtime {
    docker_url:"<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpp
    req_memory:"${mem}Gi"
  }
  output {
    File peak = "ATAC/out/~{ID}_peaks.narrowPeak"
    File peakMatrix = "ATAC/out/Peak"
    File promoterMatrix = "ATAC/out/Promoter"
    File cellReport = "ATAC/report/1.cell_report.csv"
    File mapping_csv="ATAC/report/3.mapping.csv"
    File library5_QC="ATAC/report/5.library.QC.csv"
    File plot3_DropBeadsnum="ATAC/report/plot3_DropBeadsnum.svg"
    File plot4_QC="ATAC/report/plot4_QC.svg"
    File plot7_Cluster_peak="ATAC/report/plot7_Cluster_peak.svg"
    File plot8_Cluster_depth="ATAC/report/plot8_Cluster_depth.svg"
    File plot9_Cluster_annotation="ATAC/report/plot9_Cluster_annotation.svg"
    File Table2_celltypecount="ATAC/report/Table2.celltypecount.csv"
  }
}

task cDNAanno{
        input{
        Array[String] fastq1
        Array[String] fastq2
        String barcode
        String barcode_counts_raw
        String report
        String outdir
        String refdir
        String chrSpecies
        String root
        String gtf
        Boolean? intron
        Int cpu=20
        Int mem=150
        }
        command {
                export PATH="/usr/local/bin:$PATH"
                export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib64/:$LD_LIBRARY_PATH"
                mkdir -p 02.cDNAAnno
                echo "${sep='\n' fastq1}" > 02.cDNAAnno/in1.list
                echo "${sep='\n' fastq2}" > 02.cDNAAnno/in2.list

                bcPara=02.cDNAAnno/bc.Para
                echo "in1=02.cDNAAnno/in1.list" > $bcPara
                echo "in2=02.cDNAAnno/in2.list" >> $bcPara
                echo "config=${barcode}" >> $bcPara
                echo "cbdis=02.cDNAAnno/${barcode_counts_raw}" >> $bcPara
                echo "report=02.cDNAAnno/${report}" >> $bcPara
                echo "adapter=/dsciRNA/adapter.txt" >> $bcPara
                echo "longbc" >> $bcPara
                # ${root}/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" /scRNA-seq-v3.1.5-pipeline/bin/scStar --outSAMattributes singleCell --outSAMtype BAM Unsorted --genomeDir ${refdir} --outFileNamePrefix 02.cDNAAnno/ --stParaFile $bcPara --outSAMmode NoQS --runThreadN 10

                ${root}/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/scStar_dsci  --outSAMattributes singleCell --outSAMtype BAM Unsorted --genomeDir ${refdir} --outFileNamePrefix 02.cDNAAnno/ --stParaFile $bcPara --outSAMmode NoQS --runThreadN 10 --limitOutSJcollapsed 10000000 --limitIObufferSize 350000000
                ${root}/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/Anno_dsci -I 02.cDNAAnno/Aligned.out.bam -a ${gtf} -L 02.cDNAAnno/${barcode_counts_raw} -o 02.cDNAAnno -c 10 -m chrM --species ${chrSpecies} -B ${barcode} ${true="--intron" false="" intron } --anno 1
                rm -rf 02.cDNAAnno/tmp.*.bam
                export PATH="/usr/local/bin:$PATH"
                export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib64/:/usr/local/pkgs/openssl-3.1.0-h0b41bf4_0/lib:$LD_LIBRARY_PATH"
                ${root}/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" samtools sort -T 02.cDNAAnno/tmp -@ 10 02.cDNAAnno/final.bam -o 02.cDNAAnno/final.sorted.bam
                rm -rf 02.cDNAAnno/tmp.*.bam
                export PATH="/usr/local/bin:$PATH"
                export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/pkgs/openssl-3.1.0-h0b41bf4_0/lib:/usr/local/pkgs/libtiff-4.4.0-h82bc61c_5/lib:/usr/lib64/:$LD_LIBRARY_PATH"
                ${root}/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" samtools index 02.cDNAAnno/final.sorted.bam
                mkdir -p 04.Matrix
                mkdir -p 04.Matrix/RawMatrix
                ${root}/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" /scRNA-seq-v3.1.5-pipeline/bin/PISA count -one-hit -@ 10 -cb CB -anno-tag GN -umi UB -outdir 04.Matrix/RawMatrix 02.cDNAAnno/final.sorted.bam
                rm -rf 02.cDNAAnno/Aligned.out.bam 02.cDNAAnno/final.bam

}
runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }

        output{
                File mapping = "02.cDNAAnno/alignment_report.csv"
                File annostat = "02.cDNAAnno/anno_report.csv"
                File finalsortbam = "02.cDNAAnno/final.sorted.bam"
                File beads_stat = "02.cDNAAnno/beads_stat.txt"
                File sequencing = "02.cDNAAnno/${report}"
                File rawmatrix = "04.Matrix/RawMatrix"
        }
}


task getM280UMI{
        input{
        String outdir
        String root
        String SampleID
        String finalsortbam
        String beads_stat
        String rawmatrix
        String ?ATAC_barcodeTran
        Int cpu=5
        Int mem=70
        }
        command <<<
                mkdir -p 05.MergeBead
                perl ~{root}/atcg216.pl ~{ATAC_barcodeTran} ~{root}/hex_96barcode_map_true.tsv  > 05.MergeBead/~{SampleID}_barcodeTranslate_16.txt
                awk '{print $2}' 05.MergeBead/~{SampleID}_barcodeTranslate_16.txt > 05.MergeBead/cellid.txt
                cp ~{ATAC_barcodeTran} 05.MergeBead/~{SampleID}_barcodeTranslate.txt

                awk '{print $1,$3}' ~{beads_stat} | grep -wv 0 | awk '{print $1}' > 05.MergeBead/beads_umi.txt
                awk '{print $1}' ~{ATAC_barcodeTran} > 05.MergeBead/ATAC_barcode.txt
                sort 05.MergeBead/beads_umi.txt 05.MergeBead/ATAC_barcode.txt | uniq -d > 05.MergeBead/RNA_barcodeTranslate.txt

                rm -rf 05.MergeBead/beads_umi.txt 05.MergeBead/ATAC_barcode.txt

                ~{root}/tagAdd -n 8  -bam ~{finalsortbam} -file 05.MergeBead/~{SampleID}_barcodeTranslate_16.txt -out 05.MergeBead/anno_decon.bam -tag_check CB:Z: -tag_add DB:Z: &&
                perl ~{root}/cell_stat.pl -c ~{beads_stat} -m ~{ATAC_barcodeTran} -o 05.MergeBead/ &&
                Rscript  ~{root}/CellMergeStat.R -I 05.MergeBead/~{SampleID}_barcodeTranslate.txt -O 05.MergeBead -n ~{SampleID} -b 05.MergeBead/RNA_barcodeTranslate.txt
>>>
        runtime{
             docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
             req_cpu:cpu
             req_memory:"${mem}Gi"
        }
        output{
                File anno_decon = "05.MergeBead/anno_decon.bam"
                File mergestat = "05.MergeBead/merge_cell.stat"
                File barcodeTran = "05.MergeBead/~{SampleID}_barcodeTranslate_16.txt"
                File cellid = "05.MergeBead/cellid.txt"
                File M280merge = "05.MergeBead/~{SampleID}_CellNumber_merge.png"
        }

}

task count{
  input{
    String root
    String outdir
    String anno_decon
    String ID
    String cellid
    Int cpu=5
    Int mem=80
  }
  command{
    export PATH="/usr/local/bin:$PATH"
    export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
    
    mkdir -p 04.Matrix/FilterMatrix
    /scRNA-seq-v3.1.5-pipeline/bin/PISA count -one-hit -@ 10 -cb DB -anno-tag GN -umi UB -outdir 04.Matrix/FilterMatrix ~{anno_decon} -list ~{cellid}
    /scRNA-seq-v3.1.5-pipeline/bin/PISA  attrcnt -cb DB -tags UB,GN -list ~{cellid} -dedup -o 04.Matrix/~{ID}_cell_counts.tsv -@ 20 ~{anno_decon}

  }
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpu
    req_memory:"${mem}Gi"
  }
  output{
    File matrixdir = "04.Matrix/FilterMatrix"
    File cellcount = "04.Matrix/~{ID}_cell_counts.tsv"
  }
}

task splice_matrix{
  input{
    String root
    String outdir
    String anno_decon
    String countmatrix
    String cellid
    Int cpu=5
    Int mem=100
  }
  command{
    export PATH="/usr/local/bin:$PATH"
    export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"

    mkdir -p 04.Matrix/SpliceMatrix
    mkdir -p 04.Matrix/RNAVelocityMatrix
    /scRNA-seq-v3.1.5-pipeline/bin/PISA count -one-hit -@ 10 -cb DB -ttype E,S -anno-tag GN -umi UB -outdir 04.Matrix/SpliceMatrix ~{anno_decon} -list ~{cellid}
    /scRNA-seq-v3.1.5-pipeline/bin/PISA count -one-hit -velo -@ 10 -cb DB -anno-tag GN -umi UB -outdir 04.Matrix/RNAVelocityMatrix ~{anno_decon} -list ~{cellid}

  }
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpu
    req_memory:"${mem}Gi"
  }
  output{
    File matrixdir = "04.Matrix/RNAVelocityMatrix"
    File matrixsplice= "04.Matrix/SpliceMatrix"
  }
}

task statSaturation{
  input{
    String outdir
    String ID
    String anno_decon
    String barcodeTranslate
    String beads_stat
    String root
    Int cpu=7
    Int mem=120
  }
  command{
    mkdir -p 05.MergeBead
    cp ~{barcodeTranslate} 05.MergeBead/~{ID}.barcodeMerge.tsv
    cp ~{beads_stat} 05.MergeBead/beads_stat.txt
    samtools index -c ~{anno_decon}
    python /dsciRNA/saturation_sp.py -i ~{anno_decon} -o 05.MergeBead --quality 20 --threads 10 -n ~{ID}

  }
  runtime{
    docker_url: "<registry-namespace>/<image-tag-masked-for-upload>:latest"
    req_cpu:cpu
    req_memory:"${mem}Gi"
  }
  output{
    File saturation = "05.MergeBead/~{ID}_Saturation.tsv"
    File sp = "05.MergeBead/~{ID}_Saturation.png"
  }
}