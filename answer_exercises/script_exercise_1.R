# load the required libraries
library(DSOpal)
library(dsBaseClient)
library(dsOmicsClient)

# login the connection
builder <- newDSLoginBuilder()
builder$append(server = 'study1', url = 'https://opal-demo.obiba.org', 
               user = 'dsuser', password = 'password')
builder$append(server = 'study2', url = 'https://opal-demo.obiba.org', 
               user = 'dsuser', password = 'password')

logindata <- builder$build()
conns <- datashield.login(logins = logindata)

# assign the resource to an R object
datashield.assign.resource(conns, symbol = "eSet",
                           resource = list(study1 = "workshop.GSE40732_1",
                                           study2 = "workshop.GSE40732_2"))
datashield.assign.expr(conns, symbol = "methy", expr = quote(as.resource.object(eSet)))


# answers 
ds.class("methy")
ds.nFeatures('methy')
ds.nSamples('methy')

ds.varLabels('methy')

ds.table('methy$asthma')

ans.limma <- ds.limma(model = ~ asthma,
                      Set = "methy") 
metaPvalues(ans.limma)

# logout the connection
datashield.logout(conns) 


#
# Extra ... 
#

load("c:/juan/CREAL/GitHub/brgedata/data/GSE40732_2.Rdata")
library(limma)
design <- model.matrix (~asthma, data=asthma)
fit <- lmFit(asthma, design)
fit <- eBayes(fit)
topTable(fit, coef=2)

