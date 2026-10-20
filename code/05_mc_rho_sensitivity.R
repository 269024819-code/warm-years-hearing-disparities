#!/usr/bin/env Rscript
# 补丁：ρ敏感性阶梯改用Table 1原始CI推导SE（与正文负担区间一致）
# CI-SE = (upper-lower)/3.92，来自 table1_modsevere_160.csv
library(dplyr); library(fixest)
dir.create("out", showWarnings=FALSE)
set.seed(20261008)

st <- readRDS("data/panel_stacked.rds") %>% mutate(unit=paste(iso3, female, sep="_"))
pop <- read.csv("data/country_age_sex_female_population_early_late_2003_2007_vs_2015_2019.csv")
ages <- c("55-59","60-64","65-69","70-74","75-79","80+")

# 点估计（feols，与Table 1一致）
beta_hat <- setNames(numeric(6), ages)
for (a in ages) {
  m <- feols(yld_rate ~ hot*female | iso3^female + year,
             data=st %>% filter(sev=="Moderately severe", age==a))
  beta_hat[a] <- coef(m)["hot:female"]
}

# 原始CI推导SE（正文负担区间所用）
t1 <- read.csv("data/table1_modsevere_160.csv")
se_ci <- setNames((t1$sd_hi - t1$sd_lo)/3.92, ages)
cat("CI-SE:", paste(round(se_ci,4), collapse=", "), "\n")

# 暴露与人口网格（与analysis_mc.R相同）
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
temp <- tp %>% distinct(iso3, year, .keep_all=TRUE) %>% select(iso3, year, temp_c=ts_annual_mean_temp_c)
hot <- temp %>% group_by(iso3) %>% mutate(hot=as.integer(temp_c>=quantile(temp_c,0.75,type=1))) %>% ungroup()
expo <- hot %>% mutate(period=ifelse(year<=2007,"early",ifelse(year>=2015,"late",NA))) %>%
  filter(!is.na(period)) %>% group_by(iso3, period) %>% summarise(hot_share=mean(hot), .groups="drop") %>%
  tidyr::pivot_wider(names_from=period, values_from=hot_share, names_prefix="share_")
pw <- pop %>% select(iso3, age, period, fpop=female_population_mean) %>%
  tidyr::pivot_wider(names_from=period, values_from=fpop) %>%
  rename(fpop_e=early_2003_2007, fpop_l=late_2015_2019)
grid <- pw %>% left_join(expo, by="iso3") %>% filter(!is.na(share_early))
age_idx <- match(grid$age, ages)

burden4 <- function(beta_vec) {
  b <- beta_vec[age_idx]
  c(ee=sum(b*grid$share_early*grid$fpop_e/1e5),
    le=sum(b*grid$share_late *grid$fpop_e/1e5),
    el=sum(b*grid$share_early*grid$fpop_l/1e5),
    ll=sum(b*grid$share_late *grid$fpop_l/1e5))
}
shares <- function(B) {
  se <- 0.5*((B["le"]-B["ee"])+(B["ll"]-B["el"]))
  sp <- 0.5*((B["el"]-B["ee"])+(B["ll"]-B["le"]))
  c(exposure_share=unname(se/(se+sp)), population_share=unname(sp/(se+sp)),
    B_ee=unname(B["ee"]), B_ll=unname(B["ll"]), change=unname(B["ll"]-B["ee"]))
}

N <- 20000
rho_res <- list()
for (rho in c(0, 0.3, 0.5, 0.8)) {
  Sigma <- outer(se_ci, se_ci) * (rho + (1-rho)*diag(6))
  draws <- MASS::mvrnorm(N, beta_hat, Sigma)
  ss <- t(apply(draws, 1, function(b) shares(burden4(b))))
  rho_res[[as.character(rho)]] <- data.frame(
    rho=rho,
    exposure_share_med=median(ss[,"exposure_share"]),
    exposure_share_lo=quantile(ss[,"exposure_share"],0.025),
    exposure_share_hi=quantile(ss[,"exposure_share"],0.975),
    B_ee_med=median(ss[,"B_ee"]), B_ee_lo=quantile(ss[,"B_ee"],0.025), B_ee_hi=quantile(ss[,"B_ee"],0.975),
    B_ll_med=median(ss[,"B_ll"]), B_ll_lo=quantile(ss[,"B_ll"],0.025), B_ll_hi=quantile(ss[,"B_ll"],0.975))
  cat(sprintf("rho=%.1f: share %.3f [%.3f, %.3f] | B_ee %.0f [%.0f, %.0f] | B_ll %.0f [%.0f, %.0f]\n",
      rho, median(ss[,"exposure_share"]), quantile(ss[,"exposure_share"],0.025), quantile(ss[,"exposure_share"],0.975),
      median(ss[,"B_ee"]), quantile(ss[,"B_ee"],0.025), quantile(ss[,"B_ee"],0.975),
      median(ss[,"B_ll"]), quantile(ss[,"B_ll"],0.025), quantile(ss[,"B_ll"],0.975)))
  if (rho==0) saveRDS(ss, "out/mc_draws_rho0.rds")
}
write.csv(bind_rows(rho_res), "out/table_mc_rho_sensitivity.csv", row.names=FALSE)
cat("已覆盖 out/table_mc_rho_sensitivity.csv（CI-SE版）\n")
