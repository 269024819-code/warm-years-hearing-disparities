#!/usr/bin/env Rscript
# GBD不确定性传播（修正版）：两分量模型
# 分量1（水平）：每个国家×性别×年龄序列一个系统性偏移 —— 应被FE吸收，验证β不变
# 分量2（特异）：idiosyncratic 逐年噪声，方差占比 ω ∈ {0.05, 0.10, 0.25} 做敏感性阶梯
library(dplyr); library(fixest)
dir.create("out", showWarnings=FALSE)
set.seed(20261008)
st <- readRDS("data/panel_stacked.rds") %>% mutate(unit=paste(iso3, female, sep="_"))
pop <- read.csv("data/country_age_sex_female_population_early_late_2003_2007_vs_2015_2019.csv")
ages <- c("55-59","60-64","65-69","70-74","75-79","80+")

beta_hat <- setNames(numeric(6), ages)
for (a in ages) {
  m <- feols(yld_rate ~ hot*female | iso3^female + year, data=st %>% filter(sev=="Moderately severe", age==a))
  beta_hat[a] <- coef(m)["hot:female"]
}

# 负担机器（与analysis_mc.R相同）
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
temp <- tp %>% distinct(iso3, year, .keep_all=TRUE) %>% dplyr::select(iso3, year, temp_c=ts_annual_mean_temp_c)
hot <- temp %>% group_by(iso3) %>% mutate(hot=as.integer(temp_c>=quantile(temp_c,0.75,type=1))) %>% ungroup()
expo <- hot %>% mutate(period=ifelse(year<=2007,"early",ifelse(year>=2015,"late",NA))) %>%
  filter(!is.na(period)) %>% group_by(iso3, period) %>% summarise(hot_share=mean(hot), .groups="drop") %>%
  tidyr::pivot_wider(names_from=period, values_from=hot_share, names_prefix="share_")
pw <- pop %>% dplyr::select(iso3, age, period, fpop=female_population_mean) %>%
  tidyr::pivot_wider(names_from=period, values_from=fpop) %>% rename(fpop_e=early_2003_2007, fpop_l=late_2015_2019)
grid <- pw %>% left_join(expo, by="iso3") %>% filter(!is.na(share_early))
age_idx <- match(grid$age, ages)
burden4 <- function(bv) { b <- bv[age_idx]
  c(ee=sum(b*grid$share_early*grid$fpop_e/1e5), le=sum(b*grid$share_late*grid$fpop_e/1e5),
    el=sum(b*grid$share_early*grid$fpop_l/1e5), ll=sum(b*grid$share_late*grid$fpop_l/1e5)) }
shares <- function(B) { se<-0.5*((B["le"]-B["ee"])+(B["ll"]-B["el"])); sp<-0.5*((B["el"]-B["ee"])+(B["ll"]-B["le"]))
  c(exposure_share=unname(se/(se+sp)), B_ee=unname(B["ee"]), B_ll=unname(B["ll"])) }

d_ms <- st %>% filter(sev=="Moderately severe") %>% mutate(se_gbd=(yld_upper-yld_lower)/3.92)
d_list <- lapply(ages, function(a) d_ms %>% filter(age==a))

# ── 分量1验证：水平偏移被FE吸收 ─────────────────────────────────────────────
B0 <- 200
beta_shift <- matrix(NA, B0, 6, dimnames=list(NULL, ages))
for (r in 1:B0) {
  for (j in seq_along(ages)) {
    d <- d_list[[j]]
    shift <- rnorm(length(unique(d$unit)), 0, 1)
    names(shift) <- unique(d$unit)
    d$y_pert <- pmax(0, d$yld_rate + shift[d$unit]*d$se_gbd)  # 同序列所有年份同一偏移
    beta_shift[r,j] <- coef(feols(y_pert ~ hot*female | iso3^female + year, data=d))["hot:female"]
  }
}
cat("分量1（水平偏移）: β扰动SD =", paste(round(apply(beta_shift,2,sd),4), collapse=", "),
    "（应≈0，FE吸收）\n")

# ── 分量2：idiosyncratic 敏感性阶梯 ─────────────────────────────────────────
B_draws <- 500
gbd_res <- list()
for (omega in c(0, 0.05, 0.10, 0.25)) {
  beta_g <- matrix(NA, B_draws, 6, dimnames=list(NULL, ages))
  for (r in 1:B_draws) {
    for (j in seq_along(ages)) {
      d <- d_list[[j]]
      if (omega==0) { d$y_pert <- d$yld_rate } else {
        d$y_pert <- pmax(0, d$yld_rate + rnorm(nrow(d), 0, sqrt(omega)*d$se_gbd)) }
      beta_g[r,j] <- coef(feols(y_pert ~ hot*female | iso3^female + year, data=d))["hot:female"]
    }
  }
  ss <- t(apply(beta_g, 1, function(b) shares(burden4(b))))
  gbd_res[[as.character(omega)]] <- data.frame(
    omega=omega, beta_sd_65_69=sd(beta_g[,"65-69"]),
    exposure_share_med=median(ss[,"exposure_share"]),
    exposure_share_lo=quantile(ss[,"exposure_share"],0.025, names=FALSE),
    exposure_share_hi=quantile(ss[,"exposure_share"],0.975, names=FALSE),
    B_ee_lo=quantile(ss[,"B_ee"],0.025, names=FALSE), B_ee_hi=quantile(ss[,"B_ee"],0.975, names=FALSE),
    B_ll_lo=quantile(ss[,"B_ll"],0.025, names=FALSE), B_ll_hi=quantile(ss[,"B_ll"],0.975, names=FALSE))
  cat(sprintf("omega=%.2f: β(65-69) SD=%.3f | exposure share %.3f [%.3f, %.3f] | B_ll [%.0f, %.0f]\n",
      omega, sd(beta_g[,"65-69"]), median(ss[,"exposure_share"]),
      quantile(ss[,"exposure_share"],0.025), quantile(ss[,"exposure_share"],0.975),
      quantile(ss[,"B_ll"],0.025), quantile(ss[,"B_ll"],0.975)))
  if (omega==0.10) saveRDS(list(beta=beta_g, shares=ss), "out/mc_gbd_omega10.rds")
}
write.csv(bind_rows(gbd_res), "out/table_mc_gbd_sensitivity.csv", row.names=FALSE)
cat("完成\n")
