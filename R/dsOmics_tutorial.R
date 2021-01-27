## ----load_packages------------------------------------------------------------------------------------------------------
library(DSOpal) 
library(dsBaseClient)
library(dsOmicsClient)




## ----opalOmic, echo=FALSE, fig.cap="Non-disclosive omic data analysis with DataSHIELD and Bioconductor. The figure illustrates how the `resourcer` package is used to get access to omic data through the Opal servers. Then DataSHIELD is used in the client side to perform non-disclosive data analyses.", fig.align='center'----
knitr::include_graphics("fig/dsOmics_A.jpg")


## ----omicAnal1, echo=FALSE, fig.cap="Non-disclosive omic data analysis with DataSHIELD and Bioconductor. The figure illustrates how to perform single pooled omic data analysis. The analyses are performed by using a generalized linear model (glm) on data from one or multiple sources. It makes use of `ds.glm()`, a DataSHIELD function, that uses an approach that is mathematically equivalent to placing all individual-level data from all sources in one central warehouse and analysing those data using the conventional `glm()` function in R.", fig.align='center'----
knitr::include_graphics("fig/dsOmics_B.jpg")


## ----omicAnal2, echo=FALSE, fig.cap="Non-disclosive omic data analysis with DataSHIELD and Bioconductor. The figure illustrates how to perform anlyses at genome-wide level from one or multiple sources. It runs standard Bioconductor functions at each server independently to speed up the analyses and in the case of having multiple sources, results can be meta-analyzed uning standar R functions.", fig.align='center'----
knitr::include_graphics("fig/dsOmics_C.jpg")


## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/single_study_analysis_example.Rmd", "rmd_omic/single_study_analysis_example.Rmd")


## ----child='/rmd_omic/single_study_analysis_example.Rmd', include=TRUE--------------------------------------------------

## ----pipeline_gene_expr-------------------------------------------------------------------------------------------------
builder <- newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-demo.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "RSRC.tcga_liver", driver = "OpalDriver")

logindata <- builder$build()

conns <- datashield.login(logins = logindata, assign = TRUE, 
                          symbol = "res")


## ----get_rse------------------------------------------------------------------------------------------------------------
datashield.assign.expr(conns, symbol = "rse", 
                       expr = quote(as.resource.object(res)))
ds.class("rse")


## ----dim_rse------------------------------------------------------------------------------------------------------------
ds.dim("rse")


## ----name_feature_rse---------------------------------------------------------------------------------------------------
name.features <- ds.featureNames("rse")
lapply(name.features, head)


## ----name_covar_rse-----------------------------------------------------------------------------------------------------
name.vars <- ds.featureData("rse")
lapply(name.vars, head, n=15)


## ----table_gender-------------------------------------------------------------------------------------------------------
ds.table("rse$gdc_cases.demographic.gender")


## ----voom_gender--------------------------------------------------------------------------------------------------------
ans.gender <- ds.limma(model =  ~ gdc_cases.demographic.gender, 
                       Set = "rse", type.data = "RNAseq", 
                       sva = FALSE)


## ----show_ans.gender----------------------------------------------------------------------------------------------------
ans.gender


## ----close_ds2----------------------------------------------------------------------------------------------------------
datashield.logout(conns)



## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/multiple_study_analysis_example.Rmd", "rmd_omic/multiple_study_analysis_example.Rmd")


## ----child='/rmd_omic/multiple_study_analysis_example.Rmd', include=TRUE------------------------------------------------

## ----login_assign_eSet--------------------------------------------------------------------------------------------------
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-demo.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "RSRC.GSE66351_1", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-demo.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "RSRC.GSE66351_2", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")


# Assign to the original R class (e.g ExpressionSet)
datashield.assign.expr(conns, symbol = "methy", 
                       expr = quote(as.resource.object(res)))



## ----assign_es----------------------------------------------------------------------------------------------------------
ds.class("methy")


## ----show_featureNames--------------------------------------------------------------------------------------------------
fn <- ds.featureNames("methy")
lapply(fn, head)


## ----show_phenoNames----------------------------------------------------------------------------------------------------
ds.varLabels("methy")



## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/single_cpg_analysis.Rmd", "rmd_omic/single_cpg_analysis.Rmd")


## ----child='/rmd_omic/single_cpg_analysis.Rmd', include=TRUE------------------------------------------------------------

## ----one_cpg------------------------------------------------------------------------------------------------------------
ans <- ds.lmFeature(feature = "cg07363416", 
                    model = ~ diagnosis + Sex, 
                    Set = "methy",
                    datasources = conns)
ans



## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/multiple_cpg_analysis.Rmd", "rmd_omic/multiple_cpg_analysis.Rmd")


## ----child='/rmd_omic/multiple_cpg_analysis.Rmd', include=TRUE----------------------------------------------------------

## ----multiple_cpg, eval=FALSE-------------------------------------------------------------------------------------------
## ans <- ds.lmFeature(model = ~ diagnosis + Sex,
##                     Set = "methy",
##                     datasources = conns,
##                     mc.cores = 20)


## ----limma_methy--------------------------------------------------------------------------------------------------------
ans.limma <- ds.limma(model = ~ diagnosis + Sex,
                      Set = "methy", 
                      datasources = conns)


## ----show_limma_methy---------------------------------------------------------------------------------------------------
lapply(ans.limma, head)


## ----show_annot_cols----------------------------------------------------------------------------------------------------
ds.fvarLabels("methy")


## ----remove_ans_limma, eval=FALSE, echo=FALSE---------------------------------------------------------------------------
## ds.rm("ans.limma")


## ----limma_methy_annot--------------------------------------------------------------------------------------------------
ans.limma.annot <- ds.limma(model = ~ diagnosis + Sex,
                            Set = "methy", 
                            annotCols = c("CHR", "UCSC_RefGene_Name"),
                            datasources = conns)


## ----show_limma_methy_annot---------------------------------------------------------------------------------------------
lapply(ans.limma.annot, head)


## ----meta_p-------------------------------------------------------------------------------------------------------------
ans.meta <- metaPvalues(ans.limma)
ans.meta


## ----one_cpg_val--------------------------------------------------------------------------------------------------------
res1 <- ds.lmFeature(feature = "cg13138089", 
                     model = ~ diagnosis + Sex, 
                     Set = "methy",
                     datasources = conns)
res1

res2 <- ds.lmFeature(feature = "cg13772815", 
                     model = ~ diagnosis + Sex, 
                     Set = "methy",
                     datasources = conns)
res2


## ----qqplot-------------------------------------------------------------------------------------------------------------
qqplot(ans.meta$p.meta)



## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/surrogate_variables.Rmd", "rmd_omic/surrogate_variables.Rmd")


## ----child='/rmd_omic/surrogate_variables.Rmd', include=TRUE------------------------------------------------------------

## ----login_assign_eSet_new, echo=FALSE----------------------------------------------------------------------------------
datashield.logout(conns)
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-demo.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "RSRC.GSE66351_1", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-demo.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "RSRC.GSE66351_2", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")


# Assign to the original R class (e.g ExpressionSet)
datashield.assign.expr(conns, symbol = "methy", 
                       expr = quote(as.resource.object(res)))



## ----all_cpg_sva--------------------------------------------------------------------------------------------------------
ans.sva <- ds.limma(model = ~ diagnosis + Sex, 
                    Set = "methy",
                    sva = TRUE, annotCols = c("CHR", "UCSC_RefGene_Name"))
ans.sva


## ----meta_sva-----------------------------------------------------------------------------------------------------------
ans.meta.sv <- metaPvalues(ans.sva)
ans.meta.sv


## ----close_ds-----------------------------------------------------------------------------------------------------------
datashield.logout(conns)



## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/gwas_BioC.Rmd", "rmd_omic/gwas_BioC.Rmd")


## ----child='/rmd_omic/gwas_BioC.Rmd', include=TRUE----------------------------------------------------------------------

## ----add_resources_vcf--------------------------------------------------------------------------------------------------
builder <- newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-demo.obiba.org",
               user = "dsuser", password = "password",
               resource = "RSRC.brge_vcf", driver = "OpalDriver")
logindata <- builder$build()

conns <- datashield.login(logins = logindata, assign = TRUE,
                          symbol = "res")


## ----assign_vcf---------------------------------------------------------------------------------------------------------
datashield.assign.resource(conns, symbol = "vcf.res", 
                           resource = list(study1 = "RSRC.brge_vcf"))
datashield.assign.expr(conns, symbol = "gds", 
                       expr = quote(as.resource.object(vcf.res)))


datashield.assign.resource(conns, symbol = "covars.res", 
                           resource = list(study1 = "RSRC.brge"))
datashield.assign.expr(conns, symbol = "covars", 
                       expr = quote(as.resource.data.frame(covars.res)))


## ----ls_vcf-------------------------------------------------------------------------------------------------------------
ds.ls()


## ----show_covars--------------------------------------------------------------------------------------------------------
ds.colnames("covars")


## ----show_group---------------------------------------------------------------------------------------------------------
ds.table("covars$asthma")


## ----gene_subset_1, eval=FALSE------------------------------------------------------------------------------------------
## genes <- c("A1BG","A2MP1")
## ds.getSNPSbyGen("gds", genes)


## ----gene_subset_2------------------------------------------------------------------------------------------------------
genes <- c("A1BG","A2MP1")
ds.getSNPSbyGen("gds", genes = genes, name = "subset.vcf")


## ----createGenoData-----------------------------------------------------------------------------------------------------
ds.GenotypeData(x='gds', covars = 'covars', columnId = 1, newobj.name = 'gds.Data')


## ----PCASNPs------------------------------------------------------------------------------------------------------------
ds.PCASNPS("gds", prune = TRUE)


## ----tableCovar---------------------------------------------------------------------------------------------------------
ds.table1D("covars$gender")$counts


## ----get_autosomes------------------------------------------------------------------------------------------------------
ds.getChromosomeNames("gds.Data")$autosomes


## ----HWEtest------------------------------------------------------------------------------------------------------------
ds.exactHWE("gds.Data", sexcol = "gender", male = "1", female = "2", chromosome = "22")


## ----alleleFreq---------------------------------------------------------------------------------------------------------
ds.alleleFrequency("gds.Data", sexcol = "gender", male = "1", female = "2")


## ----snp_analysis-------------------------------------------------------------------------------------------------------
ds.glmSNP(snps.fit = "rs11247693", model = asthma ~ gender + age, genoData='gds.Data')


## ----GWAS---------------------------------------------------------------------------------------------------------------
ans.bioC <- ds.GWAS('gds.Data', model=asthma~age+country)


## ----close_conns3-------------------------------------------------------------------------------------------------------
datashield.logout(conns)



## ----include = FALSE----------------------------------------------------------------------------------------------------
download.file("https://github.com/isglobal-brge/dsOmicsClient/raw/master/vignettes/rmd/plink.Rmd", "rmd_omic/plink.Rmd")


## ----child='/rmd_omic/plink.Rmd', include=TRUE--------------------------------------------------------------------------

## ----GWAS_shell_1-------------------------------------------------------------------------------------------------------
  library(DSOpal)
  library(dsBaseClient)
  library(dsOmicsClient)
  builder <- newDSLoginBuilder()
  builder$append(server = "study1", url = "https://opal-demo.obiba.org",
                 user = "dsuser", password = "password",
                 resource = "RSRC.brge_plink", driver = "OpalDriver")
  logindata <- builder$build()


## ----GWAS_shell_3-------------------------------------------------------------------------------------------------------
  conns <- datashield.login(logins = logindata, assign = TRUE,
                            symbol = "client")
  ds.class("client")


## ----GWAS_shell_4-------------------------------------------------------------------------------------------------------
  plink.arguments <- "--bfile brge --logistic --pheno brge.phe --mpheno 6 --covar brge.phe --covar-name gender,age"


## ----GWAS_shell_gwas----------------------------------------------------------------------------------------------------
  ans.plink <- ds.PLINK("client", plink.arguments)


## ----GWAS_shell_result1-------------------------------------------------------------------------------------------------
  lapply(ans.plink, names)
  
  head(ans.plink$study1$results)
  
  ans.plink$study$plink.out


## ----comparison---------------------------------------------------------------------------------------------------------
  library(tidyverse)
  # get SNP p.values (additive model - ADD)
  res.plink <- ans.plink$study1$results %>% filter(TEST=="ADD") %>%
    arrange(P)
  # compare top-10 with Biocoductor's results
  snps <- res.plink$SNP[1:10]
  plink <- res.plink %>% filter(SNP%in%snps) %>% dplyr::select(SNP, P)
  bioC <- ans.bioC$study1 %>% filter(rs%in%snps) %>% dplyr::select(rs, Score.pval)
  left_join(plink, bioC, by=c("SNP" = "rs"))


## ----maf_plink----------------------------------------------------------------------------------------------------------
  plink.arguments <- "--bfile brge --freq"
  ans.plink2 <- ds.PLINK("client", plink.arguments)
  maf.plink <- ans.plink2$study1$results
  
  plink <- maf.plink %>% filter(SNP%in%snps) %>% dplyr::select(SNP, MAF)
  bioC <- ans.bioC$study1 %>% filter(rs%in%snps) %>% dplyr::select(rs, freq)
  left_join(plink, bioC, by=c("SNP" = "rs"))


## ----close_conns4-------------------------------------------------------------------------------------------------------
  datashield.logout(conns)











## ----check--------------------------------------------------------------------------------------------------------------
library(DSOpal)
library(dsBaseClient)


builder <- newDSLoginBuilder()
builder$append(server = 'study1', url = 'https://opal-demo.obiba.org', 
               user = 'dsuser', password = 'password')
builder$append(server = 'study2', url = 'https://opal-demo.obiba.org', 
               user = 'dsuser', password = 'password')
builder$append(server = 'study3', url = 'https://opal-demo.obiba.org', 
               user = 'dsuser', password = 'password')



logindata <- builder$build()
conns <- datashield.login(logins = logindata)

datashield.assign.table(conns, symbol = "exposome",
                        table = list(study1 = "GREENSPACE.Cohort1_exposome",
                                     study2 = "GREENSPACE.Cohort2_exposome",
                                     study3 = "GREENSPACE.Cohort3_exposome"))


## -----------------------------------------------------------------------------------------------------------------------
ds.class('exposome')
ds.dim('exposome')
ds.colnames('exposome')




## ----add_resources------------------------------------------------------------------------------------------------------
datashield.assign.resource(conns, symbol = "eSet",
                           resource = list(study1 = "GREENSPACE.Cohort1_Methyl_C_blood_0y_450K",
                                           study2 = "GREENSPACE.Cohort2_Methyl_C_blood_0y_450K",
                                           study3 = "GREENSPACE.Cohort3_Methyl_C_blood_0y_450K"))

datashield.assign.expr(conns, symbol = "methy", expr = quote(as.resource.object(eSet)))


## ----check_resources----------------------------------------------------------------------------------------------------
library(dsOmicsClient)
ds.class('methy')
ds.nFeatures('methy')
ds.nSamples('methy')


## ----addtable-----------------------------------------------------------------------------------------------------------
ds.addPhenoData('methy', 'exposome', identifier = "CODE", name = 'methy2')


## ----ls-----------------------------------------------------------------------------------------------------------------
ds.ls()


## -----------------------------------------------------------------------------------------------------------------------
ds.class('methy2')


## ----showvars-----------------------------------------------------------------------------------------------------------
ds.fvarLabels('methy')
ds.fvarLabels('methy2')


## -----------------------------------------------------------------------------------------------------------------------
ds.nSamples('methy2')


## ----exwas--------------------------------------------------------------------------------------------------------------
ewas <- ds.limma(model = ~ greenyn300_0,
                      Set = "methy2")


## ----show_exwas---------------------------------------------------------------------------------------------------------
ewas


## ----meta---------------------------------------------------------------------------------------------------------------
metaPvalues(ewas)


## -----------------------------------------------------------------------------------------------------------------------
ds.fvarLabels('methy2')


## ----exwasadj-----------------------------------------------------------------------------------------------------------
ewas.adj <- ds.limma(model = ~ greenyn300_0 + sex_methyl, Set = "methy2")
metaPvalues(ewas.adj)


## ----logoutexwas--------------------------------------------------------------------------------------------------------
datashield.logout(conns)    


## -----------------------------------------------------------------------------------------------------------------------
sessionInfo()

