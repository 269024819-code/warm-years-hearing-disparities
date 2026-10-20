#!/usr/bin/env Rscript
# SE稳健性补充：two-way聚类 + wild cluster bootstrap(999) + 跨年龄合并检验 + 年龄梯度检验
library(dplyr); library(fixest)
dir.create("out", showWarnings=FALSE)
st <- readRDS("data/panel_stacked.rds") %>% mutate(unit=paste(iso3, female, sep="_"))
ages <- c("55-59","60-64","65-69","70-74","75-79","80+")
set.seed(20261008)

# ── two-way 聚类（country + year）───────────────────────────────────────────
tw <- list()
for (a in ages) {
  da <- st %>% filter(sev=="Moderately severe", age==a)
  m0 <- feols(yld_rate ~ hot*female | iso3^female + year, data=da)
  v <- vcov(m0, cluster=~iso3+year)
  tw[[a]] <- data.frame(age=a, beta=unname(coef(m0)["hot:female"]),
                        se_twoway=unname(sqrt(v["hot:female","hot:female"])))
}
tw_tab <- bind_rows(tw) %>% mutate(p_twoway=2*pnorm(-abs(beta/se_twoway)))
print(tw_tab %>% mutate(across(where(is.numeric), ~signif(.x,3))))
write.csv(tw_tab, "out/table_se_twoway.csv", row.names=FALSE)

# ── wild cluster bootstrap（Rademacher，country聚类，999次， imposing null）──
wildboot <- function(da, B=999) {
  m_u <- feols(yld_rate ~ hot*female | iso3^female + year, data=da)          # 无约束
  m_r <- feols(yld_rate ~ hot | iso3^female + year, data=da)                 # 约束（β_int=0）
  t_obs <- unname(coef(m_u)["hot:female"] / sqrt(vcov(m_u, cluster=~iso3)["hot:female","hot:female"]))
  fit_r <- fitted(m_r); res_r <- residuals(m_r)
  countries <- unique(da$iso3)
  t_boot <- numeric(B)
  for (b in 1:B) {
    w <- setNames(sample(c(-1,1), length(countries), replace=TRUE), countries)
    da$y_star <- fit_r + res_r * w[da$iso3]
    mb <- tryCatch(feols(y_star ~ hot*female | iso3^female + year, data=da), error=function(e) NULL)
    if (is.null(mb)) { t_boot[b] <- NA; next }
    vb <- tryCatch(vcov(mb, cluster=~iso3), error=function(e) NULL)
    if (is.null(vb)) { t_boot[b] <- NA; next }
    t_boot[b] <- unname(coef(mb)["hot:female"] / sqrt(vb["hot:female","hot:female"]))
  }
  t_boot <- t_boot[!is.na(t_boot) & is.finite(t_boot)]
  list(t_obs=t_obs, p=(1+sum(abs(t_boot)>=abs(t_obs)))/(length(t_boot)+1), B=length(t_boot))
}
wb_res <- list()
for (a in ages) {
  da <- st %>% filter(sev=="Moderately severe", age==a)
  r <- wildboot(da)
  wb_res[[a]] <- data.frame(age=a, t_obs=round(r$t_obs,3), p_wildboot=round(r$p,4), B=r$B)
  cat(a, "wild bootstrap p =", round(r$p,4), "\n")
}
write.csv(bind_rows(wb_res), "out/table_se_wildboot.csv", row.names=FALSE)

# ── 跨年龄合并检验（stacked all ages, age FE + country×sex FE + year FE）─────
dall <- st %>% filter(sev=="Moderately severe")
m_pool <- feols(yld_rate ~ hot*female | iso3^female + age + year, data=dall)
v_pool <- vcov(m_pool, cluster=~iso3)
b_pool <- coef(m_pool)["hot:female"]; se_pool <- sqrt(v_pool["hot:female","hot:female"])
cat(sprintf("合并检验: beta=%.3f, cluster SE=%.3f, p=%.4f\n", b_pool, se_pool, 2*pnorm(-abs(b_pool/se_pool))))

# ── 年龄梯度检验：β_a 对年龄中点的加权趋势 ──────────────────────────────────
betas <- c(0.515,0.817,1.769,2.253,2.070,4.189)
se_cl <- c(0.392,0.607,1.010,1.400,1.760,2.210)  # country聚类SE
age_mid <- c(57,62,67,72,77,82)
w <- 1/se_cl^2
gm <- lm(betas ~ age_mid, weights=w)
cat(sprintf("年龄梯度: slope=%.4f/岁, p=%.4f (加权, 聚类SE)\n", coef(gm)[2], summary(gm)$coefficients[2,4]))
# 连续交互版：hot×female×age_mid 在合并样本
dall2 <- dall %>% mutate(age_mid=recode(age,"55-59"=57,"60-64"=62,"65-69"=67,"70-74"=72,"75-79"=77,"80+"=82))
m_grad <- feols(yld_rate ~ hot*female*age_mid | iso3^female + age + year, data=dall2)
v_grad <- vcov(m_grad, cluster=~iso3)
b3 <- coef(m_grad)["hot:female:age_mid"]; s3 <- sqrt(v_grad["hot:female:age_mid","hot:female:age_mid"])
cat(sprintf("三重交互(连续年龄): beta=%.4f, cluster SE=%.4f, p=%.4f\n", b3, s3, 2*pnorm(-abs(b3/s3))))

sink("out/pooled_and_gradient_tests.txt")
cat(sprintf("Pooled Hot×Female (all ages, moderately severe): beta=%.3f, cluster SE=%.3f, p=%.4g, N=%d, countries=%d\n",
            b_pool, se_pool, 2*pnorm(-abs(b_pool/se_pool)), nrow(dall), length(unique(dall$iso3))))
cat(sprintf("Age-gradient (weighted meta-regression of 6 betas): slope=%.4f per year of age, p=%.4g\n",
            coef(gm)[2], summary(gm)$coefficients[2,4]))
cat(sprintf("Triple interaction hot×female×age_mid (continuous): beta=%.4f, cluster SE=%.4f, p=%.4g\n",
            b3, s3, 2*pnorm(-abs(b3/s3))))
sink()
cat("完成\n")
