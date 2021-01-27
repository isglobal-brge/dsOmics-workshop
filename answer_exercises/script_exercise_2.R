library(DSOpal)
library(dsBaseClient)
library(dsOmicsClient)


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
ds.class('exposome')
ds.dim('exposome')
ds.colnames('exposome')


datashield.assign.resource(conns, symbol = "eSet",
                           resource = list(study1 = "GREENSPACE.Cohort1_Methyl_C_blood_0y_450K",
                                           study2 = "GREENSPACE.Cohort2_Methyl_C_blood_0y_450K",
                                           study3 = "GREENSPACE.Cohort3_Methyl_C_blood_0y_450K"))
datashield.assign.expr(conns, symbol = "methy", expr = quote(as.resource.object(eSet)))
ds.class('methy')
ds.nFeatures('methy')
ds.nSamples('methy')

ds.addPhenoData('methy', 'exposome', identifier = "CODE", name = 'methy2')
ds.class('methy2')

ds.nSamples('methy2')


ds.fvarLabels('methy')
ds.fvarLabels('methy2')

ewas1 <- ds.limma(model = ~ greenyn300_4, Set = "methy2")
metaPvalues(ewas1)

ewas2 <- ds.limma(model = ~ ndvi100_0, Set = "methy2")
metaPvalues(ewas2)

ewas3 <- ds.limma(model = ~ ndvi100_4, Set = "methy2")
metaPvalues(ewas3)

ewas1.adj <- ds.limma(model = ~ greenyn300_4 + Bcell + CD4T + CD8T + Gran + Mono + NK, Set = "methy2")
metaPvalues(ewas1.adj)


datashield.logout(conns)    
  
  
  
  
  