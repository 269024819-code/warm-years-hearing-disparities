#!/usr/bin/env Rscript
# Fig 3 数据：国家层面三指标（暴露变化/绝对负担变化/标准化负担变化），验证计数
library(dplyr)
dir.create("out", showWarnings=FALSE)
st <- readRDS("data/panel_stacked.rds")
pop <- read.csv("data/country_age_sex_female_population_early_late_2003_2007_vs_2015_2019.csv")
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
beta <- c("55-59"=0.515,"60-64"=0.817,"65-69"=1.769,"70-74"=2.253,"75-79"=2.070,"80+"=4.189)

temp <- tp %>% distinct(iso3, year, .keep_all=TRUE) %>% dplyr::select(iso3, year, temp_c=ts_annual_mean_temp_c)
hot <- temp %>% filter(iso3 %in% unique(st$iso3)) %>% group_by(iso3) %>%
  mutate(hot=as.integer(temp_c>=quantile(temp_c,0.75,type=1))) %>% ungroup()
expo <- hot %>% mutate(period=ifelse(year<=2007,"early",ifelse(year>=2015,"late",NA))) %>%
  filter(!is.na(period)) %>% group_by(iso3, period) %>% summarise(hot_share=mean(hot), .groups="drop") %>%
  tidyr::pivot_wider(names_from=period, values_from=hot_share, names_prefix="share_") %>%
  mutate(exposure_change = share_late - share_early)

pw <- pop %>% dplyr::select(iso3, age, period, fpop=female_population_mean) %>%
  tidyr::pivot_wider(names_from=period, values_from=fpop) %>%
  rename(fpop_e=early_2003_2007, fpop_l=late_2015_2019)

g <- pw %>% left_join(expo, by="iso3") %>% mutate(beta=beta[age])
cty <- g %>% group_by(iso3) %>%
  summarise(
    exposure_change = first(exposure_change),
    B_ee = sum(beta*share_early*fpop_e/1e5),
    B_ll = sum(beta*share_late*fpop_l/1e5),
    abs_change = B_ll - B_ee,
    fpop55_e = sum(fpop_e), fpop55_l = sum(fpop_l),
    # 标准化：每10万55+女性的负担变化（用晚期人口）
    std_change_latepop = abs_change / fpop55_l * 1e5,
    # 用早晚平均人口
    std_change_meanpop = abs_change / ((fpop55_e+fpop55_l)/2) * 1e5,
    # 率的变化（不含人口增长）：β×Δshare 按晚期年龄结构加权
    rate_change = sum(beta*(share_late-share_early)*fpop_l/1e5) / fpop55_l * 1e5,
    .groups="drop")

cat("=== 计数验证 ===\n")
cat("暴露变化: 正", sum(cty$exposure_change>0), "负", sum(cty$exposure_change<0), "零", sum(cty$exposure_change==0), "(目标148/4/8)\n")
cat("绝对负担: 正", sum(cty$abs_change>0), "负", sum(cty$abs_change<0), "零", sum(cty$abs_change==0), "(目标156/4/0)\n")
cat("标准化(晚期人口): 正", sum(cty$std_change_latepop>0), "负", sum(cty$std_change_latepop<0), "(目标152/8)\n")
cat("标准化(平均人口): 正", sum(cty$std_change_meanpop>0), "负", sum(cty$std_change_meanpop<0), "\n")
cat("率变化(仅暴露): 正", sum(cty$rate_change>0), "负", sum(cty$rate_change<0), "\n")
write.csv(cty, "out/fig3_country_values.csv", row.names=FALSE)
cat("saved out/fig3_country_values.csv\n")
