---
Title: ARP2/3 complex associates with peroxisomes to participate in pexophagy in plants
Author: Jan Martinek, Petra Cifrová, Stanislav Vosolsobě, Judith García-González,
  Kateřina Malínská, Zdeňka Mauerová, Barbora Jelínková, Jana Krtková, Lenka Sikorová,
  Ian Leaves,  Imogen Sparkes, Kateřina Schwarzerová
date: "`r Sys.Date()`"
output:
  pdf_document: default
  html_document:
    df_print: paged
---


# ARP2/3 complex associates with peroxisomes to participate in pexophagy in plants

## Authors:
    
Jan Martinek, Petra Cifrová, Stanislav Vosolsobě, Judith García-González, Kateřina Malínská, Zdeňka Mauerová, Barbora Jelínková, Jana Krtková, Lenka Sikorová, Ian Leaves,  Imogen Sparkes, Kateřina Schwarzerová   
    

## Statistical analysis
  
Full R source code written by Stanislav Vosolsobě
  

### Required libraries

```{r, results='hide'}
library(betareg) # for beta regression
library(lmtest)  # for LR test after BetaReg
library(pscl)    # for zero-inflated Poisson model
library(emmeans) # for multiple comparison
```



# Figure 2b

```{r}
fig2b <- read.table("Fig2b",header = T)

m0 <- lm(data=fig2b,speed~variant)
anova(m0)

em1 <- emmeans(m0,pairwise~variant,type="response")
summary(em1)
plot(em1,type="response",comparisons = T)
```

# Figure 3e

```{r}
fig3e <- read.table("Fig3e",header = T)

m0 <- glm(data=fig3e,diffused~variant,family = quasibinomial)
anova(m0,test="LRT")

em1 <- emmeans(m0,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,type="response",comparisons = T)
```

# Figure 4a

```{r}
fig4a <- read.table("Fig4a",header = T)

m0 <- glm(data=fig4a,peroxisomes~variant,family = poisson)
anova(m0,test="LRT")

em1 <- emmeans(m0,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,comparisons = T)
```

# Figure 4b

```{r}
fig4b <- read.table("Fig4b",header = T)

m0 <- glm(data=fig4b, area~varianta,family = Gamma)
anova(m0,test="LRT")

m1 <- lm(data=fig4b, log(area)~varianta)
anova(m1)
```

# Figure 4c

```{r}
fig4c <- read.table("Fig4c",header = T)

m0 <- glm(data=fig4c, area~variant,family = Gamma)
anova(m0,test="LRT")

m1 <- lm(data=fig4c, log(area)~variant)
anova(m1)
```

# Figure 4d

```{r}
fig4d <- read.table("Fig4d",header = T)

m1 <- lm(data=fig4d, log(length) ~ variant * treatment	)
anova(m1)

em1 <- emmeans(m1, pairwise ~ treatment | variant, type="response")
summary(em1)
summary(em1)$contrasts$p.value
```

# Figure 4h

```{r}
fig4h <- read.table("Fig4h",header = T)

m0 <- glm(data=fig4h,coloc~variant,family = quasibinomial)
anova(m0,test="LRT")

em1 <- emmeans(m0,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,type="response",comparisons = T)

```

# Figure 4o

```{r}
fig4o <- read.table("Fig4o",header = T)

anova(lm(data=fig4o,tm~autophagosome))
```

# Figure S4a

```{r}
figS4a <- read.table("FigS4a",header = T)

m1 <- aov(data=figS4a,density~variant)
summary(m1)

em1 <- emmeans(m1,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,type="response",comparisons = T)

```

# Figure S4c

```{r}
figS4c <- read.table("FigS4c",header = T)

m0 <- aov(data=figS4c,fluorescence~variant)
summary(m0)

em1 <- emmeans(m0,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,type="response",comparisons = T)
```

# Figure S4d

```{r}
figS4d <- read.table("FigS4d",header = T)


m1 <- lm(data=figS4d, log(primcm)~variant*treatment)
anova(m1)

m2 <- lm(data=figS4d, log(primcm)~variant+treatment)
anova(m2)

em1 <- emmeans(m1, pairwise ~ treatment | variant, type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1, comparison=T)

```

# Figure S4e

```{r}
figS4e <- read.table("FigS4e",header = T)

m1 <- betareg(data=figS4e,I(percarea/100)~variant)
lrtest(m1)
```

# Figure 5a-d

```{r}
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

m1 <- glm(data=fig5ad,col_aut~variant,family = quasibinomial)
anova(m1,test="LRT")

```

# Figure 5e-h

```{r}
fig5eh <- read.table("Fig5efh",header = T)

m1 <- glm(data=fig5eh,round(autophagosomes)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5eh,round(arposomes)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5eh,round(colocalization)~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5eh,col_aut~variant,family = quasibinomial)
anova(m1,test="LRT")

em1 <- emmeans(m1,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,comparisons = T)

```

# Figure 5i-k

```{r}
fig5ijk <- read.table("Fig5ijk",header = T)

m1 <- glm(data=fig5ijk,arposomes~variant,family = poisson)
anova(m1,test="LRT")

m1 <- glm(data=fig5ijk,peroxisomes~variant,family = poisson)
anova(m1,test="LRT")

m1 <- betareg(data=fig5ijk,I(ap_ratio/100)~variant)
lrtest(m1)
```

# Figure S5

```{r}
figS5 <- read.table("FigS5",header = T)

m1 <- glm(data=figS5,coloc~variant,family = quasibinomial)
anova(m1,test="LRT")

em1 <- emmeans(m1,pairwise~variant,type="response")
summary(em1)
summary(em1)$contrasts$p.value
plot(em1,comparisons = T)

```
