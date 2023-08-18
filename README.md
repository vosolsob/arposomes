---
Title: ARP2/3 complex associates with peroxisomes to participate in pexophagy in plants
Author: Jan Martinek, Petra Cifrová, Stanislav Vosolsobě, Judith García-González,
  Kateřina Malínská, Zdeňka Mauerová, Barbora Jelínková, Jana Krtková, Lenka Sikorová,
  Ian Leaves,  Imogen Sparkes, Kateřina Schwarzerová
output:
  pdf_document: default
  word_document: default
---
 
# ARP2/3 complex associates with peroxisomes to participate in pexophagy in plants

### Jan Martinek, Petra Cifrová, Stanislav Vosolsobě, Judith García-González, Kateřina Malínská, Zdeňka Mauerová, Barbora Jelínková, Jana Krtková, Lenka Sikorová, Ian Leaves,  Imogen Sparkes, Kateřina Schwarzerová   
    

# Statistical analysis
  
## Full R source code written by Stanislav Vosolsobě
  

### Required libraries

```r
library(betareg) # for beta regression
library(lmtest)  # for LR test after BetaReg
library(boot)    # for bootstrapping
```

### Function for p-value estimation based on bootstrapping

```r
boot.p <- function(data,Y,X){
  prum <- function(d,i){
    g <- d[i,]
    tapply(g[,which(colnames(g)==Y)], g[,which(colnames(g)==X)], mean)
  }
  bb <- boot(data = data,prum,R=1000,
             strata = as.factor(data[,which(colnames(data)==X)]),stype = "i")
  lev <- levels(as.factor(data[,which(colnames(data)==X)]))
  A <- NULL
  B <- NULL
  pv <- rep(0,length(lev)^2-length(lev))
  meanA <- NULL
  meanB <- NULL
  for(br in 1:10){
    mn <- 1
    for(l in 1:length(lev)){
      for(m in 1:length(lev[-l])){
        lm <- (1:length(lev))[-l][m]
        if(bb$t0[lm]>bb$t0[l]){
          dci <- function(p){
            dp <- abs(bb$t0[lm] - boot.ci(bb,index=l,conf=p,type="bca")$bca[5])
            return(dp)
          }
          pval <- 1-optimize(f=dci,interval = c(0,1))$minimum
        }else{
          dci <- function(p){
            dp <- abs(bb$t0[lm] - boot.ci(bb,index=l,conf=p,type="bca")$bca[4])
            return(dp)
          }
          pval <- 1-optimize(f=dci,interval = c(0,1))$minimum
        }
        A[mn] <- lev[l]
        B[mn] <- lev[-l][m]
        pv[mn] <- pval+pv[mn]
        meanA[mn] <- bb$t0[l]
        meanB[mn] <- bb$t0[lm]
        mn <- mn + 1
      }
    }
    print(br)
  }
  out <- data.frame(A=A,B=B, pval=pv/10, meanA=meanA, meanB=meanB)
  return(out)
}
```

### Figure 2b

```r
fig2b <- read.table("Fig2b",header = T)

boot.p(fig2b,"speed","variant")

anova(lm(data=fig2b,speed~variant))
```

### Figure 3e

```r
fig3e <- read.table("Fig3e",header = T)

boot.p(fig3e,"diffused","variant")
```

### Figure 4a

```r
fig4a <- read.table("Fig4a",header = T)

boot.p(fig4a,"peroxisomes","variant")

m0 <- glm(data=fig4a,peroxisomes~variant,family = poisson)
anova(m0,test="LRT")

m1 <- glm(data=subset(fig4a,subset = fig4a$variant!="arpc2"),
          peroxisomes~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=subset(fig4a,subset = fig4a$variant!="arpc5"),
          peroxisomes~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=subset(fig4a,subset = fig4a$variant!="WT"),
          peroxisomes~variant,family = poisson)
anova(m1,test="LRT")
```

### Figure 4b

```r
fig4b <- read.table("Fig4b",header = T)

boot.p(fig4b,"area","varianta")

m0 <- glm(data=fig4b, area~varianta,family = Gamma)
anova(m0,test="LRT")

m1 <- lm(data=fig4b, log(area)~varianta)
anova(m1)
```

### Figure 4c

```r
fig4c <- read.table("Fig4c",header = T)

boot.p(fig4c,"area","variant")

m0 <- glm(data=fig4c, area~variant,family = Gamma)
anova(m0,test="LRT")

m1 <- lm(data=fig4c, log(area)~variant)
anova(m1)
```

### Figure 4d

```r
fig4d <- read.table("Fig4d",header = T)

m1 <- lm(data=fig4d, log(length) ~ variant * treatment	)
anova(m1)

vars <- as.factor(fig4d$variant)
uvars <- levels(vars)
pval <- NULL
for(i in 1:nlevels(vars)){
  dta <- fig4d[vars==uvars[i],]
  m1 <- lm(data=dta,log(length)~treatment)
  a1 <- anova(m1)
  pval[i] <- (paste(uvars[i], a1$`Pr(>F)`[1]))
  print(a1)
}
```

### Figure 4h

```r
fig4h <- read.table("Fig4h",header = T)

boot.p(fig4h,"coloc","variant")
```

### Figure 4o

```r
fig4o <- read.table("Fig4o",header = T)

boot.p(fig4o,"tm","autophagosome")

anova(lm(data=fig4o,tm~autophagosome))
```

### Figure S4a

```r
figS4a <- read.table("FigS4a",header = T)

boot.p(figS4a,"density","variant")

m1 <- aov(data=figS4a,density~variant)
summary(m1)
TukeyHSD(m1)

m1 <- lm(data=subset(figS4a,subset = figS4a$variant!="C2-K"),density~variant)
anova(m1)

m1 <- lm(data=subset(figS4a,subset = figS4a$variant!="C5-K"),density~variant)
anova(m1)

m1 <- lm(data=subset(figS4a,subset = figS4a$variant!="WT-K"),density~variant)
anova(m1)
```

### Figure S4c

```r
figS4c <- read.table("FigS4c",header = T)

boot.p(figS4c,"fluorescence","variant")

m0 <- aov(data=figS4c,fluorescence~variant)
summary(m0)
TukeyHSD(m0)

m1 <- lm(data=subset(figS4c,subset = (figS4c$variant=="WT")|
                       (figS4c$variant=="arp2")),fluorescence~variant)
anova(m1)

m1 <- lm(data=subset(figS4c,subset = (figS4c$variant=="WT")|
                       (figS4c$variant=="arpc2")),fluorescence~variant)
anova(m1)

m1 <- lm(data=subset(figS4c,subset = (figS4c$variant=="WT")|
                       (figS4c$variant=="arpc4")),fluorescence~variant)
anova(m1)

m1 <- lm(data=subset(figS4c,subset = (figS4c$variant=="WT")|
                       (figS4c$variant=="arpc5")),fluorescence~variant)
anova(m1)

m1 <- lm(data=subset(figS4c,subset = (figS4c$variant=="WT")|
                       (figS4c$variant=="atg5")),fluorescence~variant)
anova(m1)
```

### Figure S4d

```r
figS4d <- read.table("FigS4d",header = T)

boot.p(figS4d,"primcm","variant")

vars <- as.factor(unlist(lapply(strsplit(figS4d$variant,split="-"),`[`,1)))
trt <- as.factor(unlist(lapply(strsplit(figS4d$variant,split="-"),`[`,2)))
uvars <- levels(vars)
for(i in 1:length(uvars)){
  dta <- figS4d[vars==uvars[i],]
  m1 <- glm(data=dta,primcm~variant,family = poisson)
  a1 <- anova(m1,test="LRT")
  print(a1)
}

m1 <- glm(figS4d$primcm~vars*trt,family = poisson)
m2 <- glm(figS4d$primcm~vars + trt,family = poisson)
m3 <- glm(figS4d$primcm~trt,family = poisson)
anova(m1,m2,test = "LRT") 
anova(m2,m3,test = "LRT") 
anova(m3,test = "LRT") 

for(i in 1:length(uvars)){
  dta <- figS4d[vars==uvars[i],]
  m1 <- lm(data=dta,primcm~variant)
  a1 <- anova(m1)
  print(a1)
}
```

### Figure S4e

```r
figS4e <- read.table("FigS4e",header = T)

boot.p(figS4e,"percarea","variant")

m1 <- betareg(data=figS4e,I(percarea/100)~variant)
lrtest(m1)
```

### Figure 5a-d

```r
fig5ad <- read.table("Fig5ad",header = T)

# Counts were normalized by scanning areas (~ divided approx. by 2), thus the
# multiplication by 2 was necessary to prevent the loss of complexity in data
# after rounding (to get discrete values for Poisson model).
m1 <- glm(data=fig5ad,round(2*autophagosomes)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5ad,round(2*arposomes)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5ad,round(2*colocalization)~variant,family = poisson)
anova(m1,test="LRT")


boot.p(fig5ad,"autophagosomes","variant")
boot.p(fig5ad,"arposomes","variant")
boot.p(fig5ad,"colocalization","variant")
boot.p(fig5ad,"col_aut","variant")
```

### Figure 5e-h

```r
fig5eh <- read.table("Fig5efh",header = T)

m1 <- glm(data=fig5eh,round(autophagosomes)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5eh,round(arposomes)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5eh,round(colocalization)~variant,family = poisson)
anova(m1,test="LRT")


boot.p(fig5eh,"autophagosomes","variant")
boot.p(fig5eh,"arposomes","variant")
boot.p(fig5eh,"colocalization","variant")
boot.p(fig5eh,"col_aut","variant")
```

### Figure 5i-k

```r
fig5ijk <- read.table("Fig5ijk",header = T)

m1 <- glm(data=fig5ijk,arposomes~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5ijk,peroxisomes~variant,family = poisson)
anova(m1,test="LRT")

m1 <- betareg(data=fig5ijk,I(ap_ratio/100)~variant)
lrtest(m1)

boot.p(fig5ijk,"arposomes","variant")
boot.p(fig5ijk,"peroxisomes","variant")
boot.p(fig5ijk,"ap_ratio","variant")
```

### Figure S5

```r
figS5<- read.table("FigS5",header = T)

boot.p(figS5,"coloc","variant")
```
