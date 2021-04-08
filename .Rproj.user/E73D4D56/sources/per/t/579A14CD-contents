rm(list=ls())

montoya = read.csv('/Users/hhuang/Desktop/gene_exp_imm_clock/montoyaFIN_GE_coVars.csv')
kip_iage = read.csv('/Users/hhuang/Desktop/gene_exp_imm_clock/kip_iage.csv')
load('/Users/hhuang/Desktop/gene_exp_imm_clock/allExpressionData.Rd')
first_instance = read.table(file = '/Users/hhuang/Desktop/gene_exp_imm_clock/combined_data_batch_corrected.tsv', sep = '\t', header = TRUE)
first_instance[,11:ncol(first_instance)] = scale(first_instance[, 11:ncol(first_instance)])

## contains:
## * exprMats.trim- list of expression matrices from each study year.
## * pInfo        - list of donor information for each exprMat, each row corresponds to a
##                  column in exprMat.big
## * exprMat.big  - combined matrix containing all donors from all 5 years
## * pInfo.big    - combined donor info matrix, each row corresponds to a column in
##                  exprMat.big
## * exprs.qnorm  - quantile normalized version of expression matrix exprMat.big

## pathways
library(fgsea)
## kegg
kegg.pathways<- gmtPathways('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/msigDB/c2.cp.kegg.v7.1.symbols.gmt')
## GEO
geo.pathways <- gmtPathways('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/msigDB/c5.all.v7.1.symbols.gmt')


kegg_pathway_names = names(kegg.pathways)[grep('AUTOPHA', names(kegg.pathways))]
geo_pathway_names = names(geo.pathways)[grep('AUTOPHA', names(geo.pathways))]

all.kegg.genes = unlist(kegg.pathways[kegg_pathway_names])
all.geo.genes = unlist(geo.pathways[geo_pathway_names])

length(all.kegg.genes)
rel.kegg.genes = intersect(colnames(first_instance), all.kegg.genes)
rel.geo.genes = intersect(colnames(first_instance), all.geo.genes)

### linear regression ###
formula = as.formula(
  paste('immunage~', 
        paste(rel.kegg.genes, collapse = '+'),
        sep = '')
)

lfit = glm(formula = formula, data = first_instance)
### linear regression ###

library(glmnet)
library(foreach)

glmnet_coef = function(rel.gene.list){
  a <- seq(0.1, 1, 0.05)
  search <- foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(as.matrix(first_instance[,colnames(first_instance)[colnames(first_instance) %in% rel.gene.list]]),
                    first_instance[,'immunage'], 
                    family = "gaussian", nfold = 10, type.measure = "mae", paralle = TRUE, alpha = i)
    
    data.frame(cvm.1se = cv$cvm[cv$lambda == cv$lambda.1se], 
               cvm.min = cv$cvm[cv$lambda == cv$lambda.min],
               lambda.1se = cv$lambda.1se, 
               lambda.min = cv$lambda.min, 
               alpha = i)
  }
  cv3 <- search[search$cvm.min == min(search$cvm.min), ]
  md3 <- cv.glmnet(as.matrix(first_instance[,colnames(first_instance)[colnames(first_instance) %in% rel.gene.list]]),
                   first_instance[,'immunage'], 
                   family = "gaussian", nfold = 10, type.measure = "mae", paralle = TRUE, alpha = cv3$alpha)
  
  coef(md3)
  tmp_coeffs <- coef(md3, s = "lambda.1se")
  top_genes1 = data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  print(top_genes1)
  return(md3)
}

kegg.model = glmnet_coef(rel.kegg.genes)
geo.model = glmnet_coef(rel.geo.genes)

###### see where those genes in all geo belongs to
temp.selected.geo.genes = coef(geo.model, s = 'lambda.1se')
selected.geo.genes = data.frame(name = temp.selected.geo.genes@Dimnames[[1]][temp.selected.geo.genes@i + 1], coefficient = temp.selected.geo.genes@x)
selected.geo.genes = selected.geo.genes[2:nrow(selected.geo.genes),1]

all.geo.genes[which(all.geo.genes %in% selected.geo.genes)]
###### see where those genes in all geo belongs to


#####individual geo gene sets
geo.sub.models = list()
geo.sub.models <- vector(mode="list", length=length(geo_pathway_names))
names(geo.sub.models) <- geo_pathway_names
for (i in seq(1,length(geo_pathway_names))){
  print(geo_pathway_names[i])
  sub.rel.geo.genes = intersect(colnames(first_instance), unlist(geo.pathways[geo_pathway_names[i]]))
  geo.sub.models[[geo_pathway_names[[i]]]] = glmnet_coef(sub.rel.geo.genes)
}

### only these two gene sets had none zero genes
# GO_MACROAUTOPHAGY
# GO_PROCESS_UTILIZING_AUTOPHAGIC_MECHANISM
#####individual geo gene sets

##### pls #####
library(pls)
library(caret)


 
# pls.rmse.res = list()
pls.rmse.res <- vector(mode="list", length=length(geo_pathway_names))
names(pls.rmse.res) <- geo_pathway_names
for (i in seq(1:length(geo_pathway_names))){
  
  sub.rel.geo.genes.pls = intersect(colnames(first_instance), unlist(geo.pathways[geo_pathway_names[i]]))
  geo_formula = as.formula(
    paste('immunage~', 
          # paste(rel.geo.genes, collapse = '+'),
          paste(sub.rel.geo.genes.pls, collapse = '+'),
          sep = '')
  )
  
  model <- train(
    geo_formula, data = first_instance, method = "pls",
    scale = TRUE,
    trControl = trainControl("cv", number = 20),
    tuneLength = 10
  )
  # plot(model)
  
  ncomp = model$bestTune[[1]]
  
  library(plsdepot)
  # pls1 = plsreg1(first_instance[,colnames(first_instance)[colnames(first_instance) %in% rel.gene.list]],
  # pls1 = plsreg1(first_instance[,colnames(first_instance)[colnames(first_instance) %in% rel.geo.genes]],
  pls1 = plsreg1(first_instance[,colnames(first_instance)[colnames(first_instance) %in% sub.rel.geo.genes.pls]],
                 first_instance[,'immunage'],
                 comps = ncomp)
  
  # plot(first_instance$immunage, pls1$y.pred, type = 'n', xlab='Original', ylab = 'Predicted')
  # text(first_instance$immunage, pls1$y.pred, col = '#5592e3')
  rmse_res = RMSE(pls1$y.pred, first_instance$immunage)
  
  ## plotting correlation (pred vs. original)
  jpeg(paste('results/rmse/', geo_pathway_names[i], '.jpg', sep = ''))
  plot(first_instance$immunage, pls1$y.pred, xlab='Original', ylab = 'Predicted')
  title('Comparison of responses', cex.main = 0.9)
  abline(a = 0, b = 1, col = 'gray85', lwd = 2)
  text(min(first_instance$immunage) + 10, max(pls1$y.pred) - 1, paste('RMSE:', round(rmse_res, 3)))
  dev.off()
  
  # cor(first_instance$immunage, pls1$y.pred)
  
  ## correlation circle
  jpeg(paste('results/cor_circle/', geo_pathway_names[i], '.jpg', sep = ''))
  plot(pls1)
  dev.off()
  
  pls1$y.loads
  pls1$std.coefs[order(pls1$std.coefs)]
  pls1$reg.coefs[order(pls1$reg.coefs)]
  
  pls.rmse.res[[geo_pathway_names[[i]]]] = rmse_res
  
  
}
##### pls #####