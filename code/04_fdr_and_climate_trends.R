#!/usr/bin/env Rscript
# 分析C：多重比较FDR（意见3）+ 气候趋势分析（意见2）
library(dplyr); library(fixest)
dir.create("out", showWarnings=FALSE)
st <- readRDS("data/panel_stacked.rds") %>% mutate(unit=paste(iso3, female, sep="_"))
ages <- c("55-59","60-64","65-69","70-74","75-79","80+")
sevs <- c("Mild","Moderate","Moderately severe","Severe","Profound","Complete")

# ── FDR：6严重度×6年龄 Hot×Female 交互 ─────────────────────────────────────
grid <- expand.grid(sev=sevs, age=ages, stringsAsFactors=FALSE)
out <- list()
for (i in 1:nrow(grid)) {
  da <- st %>% filter(sev==grid$sev[i], age==grid$age[i])
  m <- feols(yld_rate ~ hot*female | iso3^female + year, data=da, panel.id=~unit+year)
  v <- vcov(m, "DK", lag=2)
  b <- unname(coef(m)["hot:female"]); s <- sqrt(v["hot:female","hot:female"])
  out[[i]] <- data.frame(severity=grid$sev[i], age=grid$age[i], beta=b, se_dk=s,
                         p_raw=2*pnorm(-abs(b/s)))
}
fdr <- bind_rows(out)
fdr$p_fdr <- p.adjust(fdr$p_raw, method="BH")
fdr <- fdr %>% arrange(p_raw) %>%
  mutate(sig_fdr = ifelse(p_fdr<0.05,"*",""),
         beta=round(beta,3), se_dk=round(se_dk,3), p_raw=signif(p_raw,3), p_fdr=signif(p_fdr,3))
print(as.data.frame(fdr))
write.csv(fdr, "out/table_fdr_36tests.csv", row.names=FALSE)
cat("FDR<0.05 的检验数:", sum(fdr$p_fdr<0.05), "/ 36\n")
cat("Moderately severe 各年龄 FDR:\n")
print(fdr %>% filter(severity=="Moderately severe") %>% arrange(age))

# moderately severe 是否效应最大（按年龄组内 |beta| 排名）
rk <- fdr %>% group_by(age) %>% mutate(rank_abs = rank(-abs(beta))) %>%
  filter(severity=="Moderately severe") %>% select(age, beta, rank_abs)
cat("Moderately severe 在各年龄组内 |β| 排名:\n"); print(rk)

# ── 气候趋势：160国年均温线性趋势 ───────────────────────────────────────────
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
temp <- tp %>% distinct(iso3, year, .keep_all=TRUE) %>% select(iso3, year, temp_c=ts_annual_mean_temp_c)
tr <- temp %>% group_by(iso3) %>%
  summarise(slope=coef(lm(temp_c~year))[2], p=summary(lm(temp_c~year))$coefficients[2,4], .groups="drop")
cat("\n=== 气候趋势 ===\n")
cat("160国年均温趋势: 中位 slope =", round(median(tr$slope),4), "°C/年 =",
    round(median(tr$slope)*10,3), "°C/十年\n")
cat("正趋势国家:", sum(tr$slope>0), "/", nrow(tr), "| 显著正(p<0.05):", sum(tr$slope>0 & tr$p<0.05),
    "| 显著负:", sum(tr$slope<0 & tr$p<0.05), "\n")
write.csv(tr, "out/table_climate_trends.csv", row.names=FALSE)

# Hot年份的时间分布（早期vs晚期占比，供正文）
hot <- temp %>% group_by(iso3) %>% mutate(hot=as.integer(temp_c>=quantile(temp_c,0.75,type=1))) %>% ungroup()
hd <- hot %>% summarise(early_2003_2007=mean(hot[year<=2007]), late_2015_2019=mean(hot[year>=2015]))
cat("Hot年占比: 2003-2007均值 =", round(mean(hd$early_2003_2007),3), "| 2015-2019均值 =", round(mean(hd$late_2015_2019),3), "\n")
cat("晚期Hot年占比>早期的国家:", sum(hd$late_2015_2019>hd$early_2003_2007), "/ 160\n")
