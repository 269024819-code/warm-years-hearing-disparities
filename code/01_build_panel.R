#!/usr/bin/env Rscript
# 重建分析面板：sexdiff（女−男）+ 温度 + Hot 指标，并验证 Table 1 数值
library(dplyr); library(fixest)

sg <- read.csv("data/Standard_Gender_Hearing_Panel_2003_2019.csv")
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
cmap <- read.csv("data/country_zh_en_map.csv")

sev_map <- c("轻度听力受损"="Mild","中度听力受损"="Moderate","中重度听力受损"="Moderately severe",
             "重度听力受损"="Severe","极重度听力受损"="Profound","完全听力受损"="Complete")
age_map <- c("55-59岁"="55-59","60-64岁"="60-64","65-69岁"="65-69","70-74岁"="70-74",
             "75-79岁"="75-79",">80岁"="80+")

sg2 <- sg %>%
  filter(age_group != "全部") %>%
  mutate(iso3 = cmap$iso3[match(location, cmap$zh)],
         severity_en = sev_map[severity],
         age = age_map[age_group]) %>%
  filter(!is.na(iso3))

# 女−男 差值（宽表）
w <- sg2 %>%
  select(iso3, year, age, severity_en, gender, yld_rate, yld_upper, yld_lower) %>%
  tidyr::pivot_wider(names_from = gender, values_from = c(yld_rate, yld_upper, yld_lower)) %>%
  mutate(sexdiff = yld_rate_女 - yld_rate_男,
         # 差值的GBD不确定性（独立性假设下合成）
         sexdiff_se = sqrt(((yld_upper_女 - yld_lower_女)/3.92)^2 + ((yld_upper_男 - yld_lower_男)/3.92)^2))

# 温度（国家-年层面，167国）+ Hot 指标
temp <- tp %>% distinct(iso3, year, .keep_all = TRUE) %>%
  select(iso3, year, temp_c = ts_annual_mean_temp_c, anomaly = ts_temp_anomaly_vs_2003_2012_c)
hot <- temp %>% group_by(iso3) %>%
  mutate(hot = as.integer(temp_c >= quantile(temp_c, 0.75, type = 1))) %>% ungroup()

panel <- w %>% inner_join(hot, by = c("iso3","year"))
cat("合并后面板:", nrow(panel), "行 | 国家数:", length(unique(panel$iso3)), "\n")
cat("每国Hot年数:", panel %>% filter(severity_en=="Moderately severe", age=="65-69") %>%
      group_by(iso3) %>% summarise(h=sum(hot)) %>% pull(h) %>% unique(), "\n")

saveRDS(panel, "data/panel_sexdiff.rds")

# ── 验证：moderately severe × 65-69 回归，β 应 ≈ 1.769（Table 1 记录值）────────
d <- panel %>% filter(severity_en == "Moderately severe")
for (a in c("55-59","60-64","65-69","70-74","75-79","80+")) {
  da <- d %>% filter(age == a)
  m <- feols(sexdiff ~ hot | iso3 + year, data = da)
  cat(sprintf("%s: beta=%.3f (n=%d, countries=%d)\n", a, coef(m)["hot"], nrow(da), length(unique(da$iso3))))
}
