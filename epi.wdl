task OptitypeDNA {

    File read1
    File read2

    String baseName = basename(read1, "_1P.fastq.gz")

    command {
	    python /usr/local/bin/OptiType/OptiTypePipeline.py --dna -i ${read1} ${read2} -o $(pwd)/${baseName}.hla.tsv;
        cat `find . -name '*result.tsv' -print`
    }

    runtime {
	    docker: "fred2/optitype:release-v1.3.1"
    }

    output {
	    File hlaTypes = stdout()
    }

}

task Kallisto {

    File index
    Int fragmentLength
    Int standardDev
    Int nThreads
    File rnaInput1
    File rnaInput2

    command {
        kallisto quant -i ${index} -l ${fragmentLength} -s ${standardDev} -o . -t ${nThreads} \
        ${rnaInput1} ${rnaInput2}
    }

    runtime {
        cpu: "${nThreads}"
        docker: "insilicodb/kallisto:1.0.0"
    }

    output {
        File out = "abundance.tsv"
    }

}

task TrimmomaticPE {
    Int nThreads
    File read1
    File read2
    File adapters
    Int seedMismatches
    Int palindromeClipThreshold
    Int simpleClipThreshold
    Int crop
    Int leading
    Int trailing
    Int windowSize
    Int requiredQuality
    Int minLen
    String baseName


    runtime {
        cpu: "${nThreads}"
        docker: "comics/trimmomatic:0.36"
    }

    command {
        java -jar $TRIMMOMATIC PE -trimlog ${baseName}.trimlog.txt -threads ${nThreads} \
        ${read1} ${read2} \
        -baseout ${baseName}.fastq.gz \
        ILLUMINACLIP:${adapters}:${seedMismatches}:${palindromeClipThreshold}:${simpleClipThreshold} CROP:${crop} \
        LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${windowSize}:${requiredQuality} MINLEN:${minLen}
    }

    output {
        File out1Paired = "${baseName}_1P.fastq.gz"
        File out1Unpaired = "${baseName}_1U.fastq.gz"
        File out2Paired = "${baseName}_2P.fastq.gz"
        File out2Unpaired = "${baseName}_2U.fastq.gz"
        File trimlog = "${baseName}.trimlog.txt"
    }
}

task BWAMem {
    File read1
    File read2
    File refGenomeFolder
    String refGenome
    Int nThreads

    String baseName = basename(read1, "_1P.fastq.gz")

    command {
        bwa mem -M -t ${nThreads} ${refGenomeFolder}/${refGenome} ${read1} ${read2} > ${baseName}.bwa.sam
    }

    runtime {
        cpu: "${nThreads}"
        docker: "biocontainers/bwa:v0.7.17-3-deb_cv1"
    }

    output {
        File outputAligned = "${baseName}.bwa.sam"
    }
}

task STAR {

    File rnaFastqGz1
    File rnaFastqGz2
    File genomeDir

    String baseName = basename(rnaFastqGz1, "_reads1.fastq.gz")

    Int numThreads

    command {
        STAR --runMode alignReads --runThreadN ${numThreads} --genomeDir ${genomeDir} \
        --readFilesIn ${rnaFastqGz1} ${rnaFastqGz2} --readFilesCommand zcat --outFileNamePrefix ${baseName}. \
        --outSAMtype BAM Unsorted --outReadsUnmapped Fastx
    }

    runtime {
        cpu: "${numThreads}"
        docker: "dceoy/star:latest"
    }

    output {
        File rnaAligned = "${baseName}.Aligned.out.bam"
    }
}

task Samtools {
    File alignedSam
    String baseName

    command {
        samtools sort ${alignedSam} > ${baseName}.sorted.bam;
        samtools index ${baseName}.sorted.bam > ${baseName}.sorted.bam.bai
    }

    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    output {
        File outputSorted = "${baseName}.sorted.bam"
        File samIndex = "${baseName}.sorted.bam.bai"
    }

}

task Picard {

    File sortedBam
    String readGroupName

    String baseName = basename(sortedBam, ".sorted.bam")

    command {
        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /usr/picard/picard.jar  \
            MarkDuplicates I=${sortedBam} O=${baseName}.bam M=${baseName}.txt;
        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /usr/picard/picard.jar  \
            AddOrReplaceReadGroups I=${baseName}.bam O=${baseName}-RG.bam \
            RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${readGroupName};
        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2  -jar /usr/picard/picard.jar \
            BuildBamIndex I=${baseName}-RG.bam O=${baseName}-RG.bam.bai
    }

    runtime {
        cpu: 4
        docker: "broadinstitute/picard:2.21.2"
    }

    output {
        File bamReadGroup = "${baseName}-RG.bam"
        File bamReadGroupIndex = "${baseName}-RG.bam.bai"
    }
}

task BQSR {

    File inputBamRG
    File referenceGenomeFolder
    String referenceGenomeName
    File knownSitesFolder
    String snp
    String goldIndels
    String knownIndels

    String baseName = basename(inputBamRG, "-RG.bam")

    command {
        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /gatk/gatk.jar BaseRecalibrator \
            -R ${referenceGenomeFolder}/${referenceGenomeName} -I ${inputBamRG} \
            --known-sites ${knownSitesFolder}/${snp} --known-sites ${knownSitesFolder}/${goldIndels} \
            --known-sites ${knownSitesFolder}/${knownIndels} -O ${baseName}_recall_data.table;

        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /gatk/gatk.jar ApplyBQSR \
            -reference ${referenceGenomeFolder}/${referenceGenomeName} -input ${inputBamRG} \
            -bqsr ${baseName}_recall_data.table -output ${baseName}-BQSR.bam;

        mv ${baseName}-BQSR.bai ${baseName}-BQSR.bam.bai;
    }

    runtime {
        cpu: 4
        docker: "broadinstitute/gatk:4.1.4.0"
    }

    output {
        File recallTable = "${baseName}_recall_data.table"
        File bamRecalibrated = "${baseName}-BQSR.bam"
        File bamRecalibratedIndex = "${baseName}-BQSR.bam.bai"
    }
}

task HaplotypeCaller {

    File bamRecalibrated
    File bamRecalibratedIndex
    File refGenomeFolder
    String refGenome

    String baseName = basename(bamRecalibrated, "-N-WEX-BQSR.bam")

    command {
        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /gatk/gatk.jar \
            HaplotypeCaller -R ${refGenomeFolder}/${refGenome} -I ${bamRecalibrated} \
            -O ${baseName}.germline.vcf
    }

    runtime {
        cpu: 4
        docker: "broadinstitute/gatk:4.1.4.0"
    }

    output {
        File variants = "${baseName}.germline.vcf"
    }

}

task Strelka {

    File normalBam
    File normalBamIndex
    File tumorBam
    File tumorBamIndex
    File refGenomeFolder
    String refGenome
    Int nCores

    command <<<
	    /opt/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam=${normalBam} --tumorBam=${tumorBam} \
            --referenceFasta=${refGenomeFolder}/${refGenome} --exome --runDir=$(pwd);

        python runWorkflow.py -m local -j ${nCores};
        zcat results/variants/somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > results/variants/somatic.snvs.gt.vcf;
        zcat results/variants/somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > results/variants/somatic.indels.gt.vcf;
    >>>

    runtime {
	    cpu: "${nCores}"
	    docker: "mgibio/strelka:2.9.9"
    }

    output {
	    File variantsSNP = "results/variants/somatic.snvs.gt.vcf"
	    File variantsIndel = "results/variants/somatic.indels.gt.vcf"
    }
}

task Mutect2 {

    File tumorRecalibrated
    File tumorRecalibratedIndex
    File normalRecalibrated
    File normalRecalibratedIndex
    File refGenomeFolder
    String refGenome

    String baseName = basename(normalRecalibrated, "-N-WEX-BQSR.bam")

    command {
        java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /gatk/gatk.jar Mutect2 \
            -R ${refGenomeFolder}/${refGenome} -I ${tumorRecalibrated} -I ${normalRecalibrated} \
            --normal Normal --independent-mates --native-pair-hmm-threads 8 -O mutect2.vcf.gz
    }

    runtime {
        cpu: 4
        docker: "broadinstitute/gatk:4.1.4.0"
    }

    output {
        File variants = "mutect2.vcf.gz"
    }

}

task SamtoolsMpileUp {

    File normalBam
    File tumorBam
    File refGenomeFolder
    String refGenome

    String baseName = basename(normalBam, "-BQSR.bam")

    command {
    	samtools mpileup -f ${refGenomeFolder}/${refGenome} ${normalBam} ${tumorBam} -o ${baseName}.mpileup
    }

    runtime {
	    docker: "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    output {
	    File pileUp = "${baseName}.mpileup"
    }
}

task Varscan {

    File inputMpileUp
    String baseName = basename(inputMpileUp, ".mpileup")

    command {
	    java -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -jar /opt/varscan/VarScan.jar somatic \
            ${inputMpileUp} ${baseName} --mpileup 1 --output-vcf 1;
    }

    runtime {
	    cpu: 4
	    docker: "mgibio/varscan:v2.4.2"
    }

    output {
	    File variantsSNP = "${baseName}.snp.vcf"
	    File variantsIndel = "${baseName}.indel.vcf"
    }
}

task BcfToolsConcat {

    File snp
    File indel
    String outName

    command {
        bgzip -f ${snp};
        bgzip -f ${indel};

        tabix -f ${snp}.gz;
        tabix -f ${indel}.gz;

        bcftools concat -a -o ${outName}.vcf.gz -O z ${snp}.gz ${indel}.gz;
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    }

    output {
        File vcfConcat = "${outName}.vcf.gz"
    }
}

task BcfToolsNorm {

    File vcfFile
    File refGenomeFolder
    String refGenome

    String baseName = basename(vcfFile, ".vcf.gz")

    command {
        bcftools sort ${vcfFile} -O z -o ${baseName}.sorted.vcf.gz;
        bcftools norm --check-ref wx -f ${refGenomeFolder}/${refGenome} -m - ${baseName}.sorted.vcf.gz -O z -o ${baseName}.norm.vcf.gz
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    }

    output {
        File vcfNorm = "${baseName}.norm.vcf.gz"
    }
}

task Epi {

    File dna
    File dnaIndex
    File rna
    File rnaIndex
    File vcfMutect
    File vcfStrelka
    File vcfVarscan
    File vcfGermline
    File proteins
    File gtf
    File refGenomeFolder
    String refGenome
    File hla
    File tpm
    File vepCacheDir

    command {
        Epi --dna ${dna} --rna ${rna} --somatic ${vcfMutect} ${vcfStrelka} ${vcfVarscan} \
        --vc "Mutect2" "Strelka" "Varscan" --germline ${vcfGermline} --proteins ${proteins} \
        --gtf ${gtf} --ref ${refGenomeFolder}/${refGenome} --hla ${hla} --tpm ${tpm} \
        --phaser "/opt/phaser/phaser" --vepCacheDir ${vepCacheDir}
    }

    runtime {
        cpu: 4
        docker: "epi"
    }

    output {
        File results = "epitopes.tsv"
    }

}


workflow EpitopePrediction {

    String seqNormal1
    String seqNormal2
    String seqTumor1
    String seqTumor2
    String seqRNA1
    String seqRNA2

    String seqAdapters
    Int cropLength

    String refGenomeFolder
    String refGenomeName
    String refGenomeBWAFolder
    String refGenomeBWAPrefix

    String knownSitesFolder
    String refSNP
    String refGoldIndels
    String refKnownIndels
    String cdnaIndex
    String proteins
    String gtf
    String STARGenomeDir
    String vepCacheDir

    # TPM count #
    call Kallisto as KallistoRNA {
        input:
            index = cdnaIndex,
            fragmentLength = 180,
            standardDev = 20,
            nThreads  = 10,
            rnaInput1 = seqRNA1,
            rnaInput2 = seqRNA2
    }


    # Preprocessing #
    call TrimmomaticPE as TrimmoNormal {
        input:
            nThreads = 4,
            read1 = seqNormal1,
            read2 = seqNormal2,
            adapters = seqAdapters,
            seedMismatches = 2,
            palindromeClipThreshold = 30,
            simpleClipThreshold = 10,
            crop = cropLength,
            leading = 3,
            trailing = 3,
            windowSize = 4,
            requiredQuality = 15,
            minLen = 36,
            baseName = "DNA-N"
    }

    call TrimmomaticPE as TrimmoTumor {
        input:
            nThreads = 4,
            read1 = seqTumor1,
            read2 = seqTumor2,
            adapters = seqAdapters,
            seedMismatches = 2,
            palindromeClipThreshold = 30,
            simpleClipThreshold = 10,
            crop = cropLength,
            leading = 3,
            trailing = 3,
            windowSize = 4,
            requiredQuality = 15,
            minLen = 36,
            baseName = "DNA-T"
    }

    # HLA Typing #
    call OptitypeDNA {
        input:
            read1 = TrimmoNormal.out1Paired,
            read2 = TrimmoNormal.out2Paired,
    }

    # Alignment #
    call BWAMem as BWAMemNormal {
        input:
            nThreads = 20,
            read1 = TrimmoNormal.out1Paired,
            read2 = TrimmoNormal.out2Paired,
            refGenomeFolder = refGenomeBWAFolder,
            refGenome = refGenomeBWAPrefix
    }

    call BWAMem as BWAMemTumor {
        input:
            nThreads = 20,
            read1 = TrimmoTumor.out1Paired,
            read2 = TrimmoTumor.out2Paired,
            refGenomeFolder = refGenomeBWAFolder,
            refGenome = refGenomeBWAPrefix
    }

    call STAR {
        input:
            rnaFastqGz1 = seqRNA1,
            rnaFastqGz2 = seqRNA2,
            genomeDir = STARGenomeDir,
            numThreads = 10
    }


    # Sort and Index #
    call Samtools as SamToolsNormal {
        input:
            alignedSam = BWAMemNormal.outputAligned,
            baseName = "DNA-N"
    }

    call Samtools as SamToolsTumor {
        input:
            alignedSam = BWAMemTumor.outputAligned,
            baseName = "DNA-T"
    }

    call Samtools as SamToolsRNA {
        input:
            alignedSam = STAR.rnaAligned,
            baseName = "RNA-T"
    }


    # Duplicate Marking #
    call Picard as PicardNormal {
        input:
            sortedBam = SamToolsNormal.outputSorted,
            readGroupName = "Normal"
    }

    call Picard as PicardTumor {
        input:
            sortedBam = SamToolsTumor.outputSorted,
            readGroupName = "Tumor"
    }

    call Picard as PicardRNA {
        input:
            sortedBam = SamToolsRNA.outputSorted,
            readGroupName = "Tumor"
    }


    # Base Quality Score Recalibration #
    call BQSR as BQSRNormal {
        input:
            inputBamRG = PicardNormal.bamReadGroup,
            referenceGenomeFolder = refGenomeFolder,
            referenceGenomeName = refGenomeName,
            knownSitesFolder = knownSitesFolder,
            snp = refSNP,
            goldIndels = refGoldIndels,
            knownIndels = refKnownIndels
    }

    call BQSR as BQSRTumor {
        input:
            inputBamRG = PicardTumor.bamReadGroup,
            referenceGenomeFolder = refGenomeFolder,
            referenceGenomeName = refGenomeName,
            knownSitesFolder = knownSitesFolder,
            snp = refSNP,
            goldIndels = refGoldIndels,
            knownIndels = refKnownIndels
    }


    # Variant Calling #
    call HaplotypeCaller {
        input:
            bamRecalibrated = BQSRNormal.bamRecalibrated,
            bamRecalibratedIndex = BQSRNormal.bamRecalibratedIndex,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }

    call Strelka {
        input:
            normalBam = BQSRNormal.bamRecalibrated,
            normalBamIndex = BQSRNormal.bamRecalibratedIndex,
            tumorBam = BQSRTumor.bamRecalibrated,
            tumorBamIndex = BQSRTumor.bamRecalibratedIndex,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName,
            nCores = 8
    }

    call Mutect2 {
        input:
            normalRecalibrated = BQSRNormal.bamRecalibrated,
            normalRecalibratedIndex = BQSRNormal.bamRecalibratedIndex,
            tumorRecalibrated = BQSRTumor.bamRecalibrated,
            tumorRecalibratedIndex = BQSRTumor.bamRecalibratedIndex,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }

    call SamtoolsMpileUp {
        input:
            normalBam = BQSRNormal.bamRecalibrated,
            tumorBam = BQSRTumor.bamRecalibrated,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }

    call Varscan {
        input:
            inputMpileUp = SamtoolsMpileUp.pileUp
    }


    # Postprocess Variant Caller Output #
    call BcfToolsConcat as BcfToolsConcatStrelka {
        input:
            snp = Strelka.variantsSNP,
            indel = Strelka.variantsIndel,
            outName = "strelka"
    }

    call BcfToolsConcat as BcfToolsConcatVarscan {
        input:
            snp = Varscan.variantsSNP,
            indel = Varscan.variantsIndel,
            outName = "varscan"
    }

    call BcfToolsNorm as BcfToolsNormStrelka {
        input:
            vcfFile = BcfToolsConcatStrelka.vcfConcat,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }

    call BcfToolsNorm as BcfToolsNormMutect2 {
        input:
            vcfFile = Mutect2.variants,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }

    call BcfToolsNorm as BcfToolsNormVarscan {
        input:
            vcfFile = BcfToolsConcatVarscan.vcfConcat,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }

    call BcfToolsNorm as BcfToolsNormHaplotypeCaller {
        input:
            vcfFile = HaplotypeCaller.variants,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName
    }


    # Neoeptiope Prediction #
    call Epi {
        input:
            dna = BQSRTumor.bamRecalibrated,
            dnaIndex = BQSRTumor.bamRecalibratedIndex,
            rna = PicardRNA.bamReadGroup,
            rnaIndex = PicardRNA.bamReadGroupIndex,
            vcfMutect = BcfToolsNormMutect2.vcfNorm,
            vcfStrelka = BcfToolsNormStrelka.vcfNorm,
            vcfVarscan = BcfToolsNormVarscan.vcfNorm,
            vcfGermline = BcfToolsNormHaplotypeCaller.vcfNorm,
            proteins = proteins,
            gtf = gtf,
            refGenomeFolder = refGenomeFolder,
            refGenome = refGenomeName,
            hla = OptitypeDNA.hlaTypes,
            tpm = KallistoRNA.out,
            vepCacheDir = vepCacheDir
    }
}
