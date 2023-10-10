nextflow.enable.dsl=2

process create_sparse_grm {
    tag "Build GRM on ${dataset}"
    //publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(bed), path(bim), path(fam)
    output:
       tuple val(dataset), 
            file("${dataset}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"),
               file("${dataset}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
    script:
        base = bed.baseName

    """
    Rscript /usr/local/bin/createSparseGRM.R  \
    --plinkFile=${base} \
    --nThreads=72  \
    --outputPrefix=${dataset}       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05
    """
}

process saige_logistic_step_1 {
    tag "Run logistic regression on ${dataset}_${subpop}"
    publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(bed), path(bim), path(fam), val(pheno_label), path(phenotype), path(sparseGRM_mtx), path(sparseGRM_mtx_sample)

    output:
       tuple val(dataset), val(pheno_label), 
            file("${dataset}_${pheno_label}_saige_out.rda"),
               file("${dataset}_${pheno_label}_saige_out.varianceRatio.txt")

    script:
        out = "${dataset}_${pheno_label}_saige_out"
        base = bed.baseName

        """
        Rscript /usr/local/bin/step1_fitNULLGLMM.R \\
        --sparseGRMFile=${sparseGRM_mtx} \\
        --sparseGRMSampleIDFile=${sparseGRM_mtx_sample} \\
        --useSparseGRMtoFitNULL=TRUE    \\
        --plinkFile=${base} \\
        --phenoFile=${phenotype} \\
        --skipVarianceRatioEstimation=FALSE \\
        --phenoCol=${pheno_label}_pheno \\
        --covarColList=sex,age,age_squared,age_sex \\
        --sampleIDColinphenoFile=IID \\
        --traitType=binary \\
        --outputPrefix=${out} \\
        --nThreads=72 \\
        --IsOverwriteVarianceRatioFile=TRUE
        """
}

process saige_logistic_step_2 {
    tag "Run logistic regression on ${dataset}_${subpop}"
    publishDir "${params.outdir}/Regression_results/SAIGE", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(vcf), path(index), path(id), val(pheno_label), path(rda), path(varianceRatio), path(sparseGRM_mtx), path(sparseGRM_mtx_sample)

    output:
       tuple val(dataset), val(pheno_label), file(out)

    script:
        out = "${dataset}_${pheno_label}_binary.SAIGE.vcf.genotype.txt"

        """
        Rscript /usr/local/bin/step2_SPAtests.R \\
        --vcfFile=${vcf} \\
        --vcfFileIndex=${index} \\
        --vcfField=GT --chrom=6 \\
        --minMAF=0.01 --minMAC=1 \\
        --sampleFile=${id} \\
        --GMMATmodelFile=${rda} \\
        --varianceRatioFile=${varianceRatio} \\
        --SAIGEOutputFile=${out} \\
        --LOCO=FALSE \\
        --AlleleOrder=ref-first \\
        --is_Firth_beta=TRUE  \\
        --pCutoffforFirth=0.1 \\
        --sparseGRMFile=${sparseGRM_mtx} \\
        --sparseGRMSampleIDFile=${sparseGRM_mtx_sample} \\
        --is_fastTest=TRUE
        """
}

workflow{
     //step 0
     plink_ch = Channel.fromList(params.whole_plink)
     create_sparse_grm(plink_ch)

     //step 1
     //cov_pheno_ch = Channel.fromList(params.cov_pheno)
     //plink_ch_1 = Channel.fromList(params.whole_plink)
     //input = plink_ch_1
       // .combine(cov_pheno_ch)
        //.combine(create_sparse_grm.out)
     //saige_logistic_step_1(input)
  
     //step2
    //step1_out = saige_logistic_step_1.out
    //vcf_ch = Channel.fromList(params.chr6_vcf)
    //ids_ch = Channel.fromPath(params.whole_ids)
    //reg_2_ch = vcf_ch.combine(ids_ch)
      //    .combine(step1_out, by:0)
    //saige_logistic_step_2(reg_2_ch)

}
