#!/usr/bin/env Rscript
# SE稳健性（老师意见7）：DK(2) / country聚类 / region×year FE / leave-one-region-out
library(dplyr); library(fixest); library(countrycode)
dir.create("out", showWarnings = FALSE)

st <- readRDS("data/panel_stacked.rds")
regions <- data.frame(iso3=unique(st$iso3), region=countrycode(unique(st$iso3), "iso3c", "region"))
st <- st %>% left_join(regions, by="iso3") %>% mutate(unit = paste(iso3, female, sep="_"))
cat("区域映射缺失:", sum(is.na(st$region)), "| 区域数:", length(unique(st$region)), "\n")

ages <- c("55-59","60-64","65-69","70-74","75-79","80+")
res <- list(); loro <- list()
for (a in ages) {
  da <- st %>% filter(sev=="Moderately severe", age==a)
  m0 <- feols(yld_rate ~ hot*female | iso3^female + year, data=da, panel.id=~unit+year)
  dk <- vcov(m0, "DK", lag=2)
  cl <- vcov(m0, cluster=~iso3)
  m_ry <- feols(yld_rate ~ hot*female | iso3^female + region^year, data=da)
  res[[a]] <- data.frame(age=a, beta=unname(coef(m0)["hot:female"]),
    se_dk=sqrt(dk["hot:female","hot:female"]),
    se_cluster=unname(sqrt(cl["hot:female","hot:female"])),
    beta_regionyear=unname(coef(m_ry)["hot:female"]),
    se_regionyear=unname(sqrt(vcov(m_ry, cluster=~iso3)["hot:female","hot:female"])))
  for (rg in na.omit(unique(da$region))) {
    m <- feols(yld_rate ~ hot*female | iso3^female + year, data=da %>% filter(region!=rg))
    loro[[paste(a,rg)]] <- data.frame(age=a, dropped_region=rg, beta=unname(coef(m)["hot:female"]),
                                      n_countries=length(unique(da$iso3[da$region!=rg])))
  }
}
se_tab <- bind_rows(res) %>% mutate(
  p_dk=2*pnorm(-abs(beta/se_dk)), p_cluster=2*pnorm(-abs(beta/se_cluster)),
  p_regionyear=2*pnorm(-abs(beta_regionyear/se_regionyear)))
print(se_tab %>% mutate(across(where(is.numeric), ~signif(.x,3))))
write.csv(se_tab, "out/table_se_robustness.csv", row.names=FALSE)

loro_tab <- bind_rows(loro)
write.csv(loro_tab, "out/table_leave_one_region_out.csv", row.names=FALSE)
print(loro_tab %>% group_by(age) %>% summarise(min_beta=round(min(beta),3), max_beta=round(max(beta),3), .groups="drop"))
cat("SE稳健性完成\n")
