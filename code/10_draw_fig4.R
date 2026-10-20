#!/usr/bin/env Rscript
# Fig 4 v9（老师意见）：左=waterfall（保留），右=水平贡献条+95%区间，"modeled contribution share"
library(ggplot2); library(dplyr); library(patchwork); library(MASS)
dir.create("out", showWarnings=FALSE)
set.seed(20261008)

theme_nature <- function(base_size=9) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
          panel.grid.major.x=element_line(color="#E2E8F0", linewidth=0.4),
          axis.line=element_line(color="#4A5568", linewidth=0.4),
          axis.ticks=element_line(color="#4A5568", linewidth=0.3),
          axis.title=element_text(color="#2D3748", size=9),
          axis.text=element_text(color="#4A5568", size=8),
          plot.tag=element_text(face="bold", size=12),
          text=element_text(family="Liberation Sans"))
}

# ── 重建负担机器 + 完整MC draws（含四角与Shapley分量）──────────────────────
st <- readRDS("data/panel_stacked.rds")
pop <- read.csv("data/country_age_sex_female_population_early_late_2003_2007_vs_2015_2019.csv")
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
ages <- c("55-59","60-64","65-69","70-74","75-79","80+")
beta_hat <- c(0.515,0.817,1.769,2.253,2.070,4.189)
se_dk <- c(0.222,0.232,0.4447,0.6314,0.8834,1.1515)  # 原始报告CI反推SE

temp <- tp %>% distinct(iso3, year, .keep_all=TRUE) %>% dplyr::select(iso3, year, temp_c=ts_annual_mean_temp_c)
hot <- temp %>% filter(iso3 %in% unique(st$iso3)) %>% group_by(iso3) %>%
  mutate(hot=as.integer(temp_c>=quantile(temp_c,0.75,type=1))) %>% ungroup()
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

N <- 20000
draws <- MASS::mvrnorm(N, beta_hat, diag(se_dk^2))  # ρ=0 主情形
mc <- t(apply(draws, 1, function(b) {
  B <- burden4(b)
  se <- 0.5*((B["le"]-B["ee"])+(B["ll"]-B["el"])); sp <- 0.5*((B["el"]-B["ee"])+(B["ll"]-B["le"]))
  c(B_ee=unname(B["ee"]), B_el=unname(B["el"]), B_le=unname(B["le"]), B_ll=unname(B["ll"]),
    sh_exp=unname(se), sh_pop=unname(sp),
    exp_share=unname(se/(se+sp)), pop_share=unname(sp/(se+sp)))
}))
saveRDS(mc, "out/mc_draws_full.rds")
q <- function(x) round(quantile(x, c(0.025,0.5,0.975)), 3)
cat("B_ee:", q(mc[,"B_ee"]), "| B_ll:", q(mc[,"B_ll"]), "\n")
cat("sh_exp:", q(mc[,"sh_exp"]), "| sh_pop:", q(mc[,"sh_pop"]), "\n")
cat("exp_share:", q(mc[,"exp_share"]), "| pop_share:", q(mc[,"pop_share"]), "\n")

# ── 面板 a：waterfall ───────────────────────────────────────────────────────
pt <- apply(mc, 2, median)
wf <- data.frame(
  x = 1:4,
  label = c("Early period\n(2003-2007)", "+ Warm-year\nexposure", "+ Female\npopulation ≥55", "Late period\n(2015-2019)"),
  ymin = c(0, pt["B_ee"], pt["B_ee"]+pt["sh_exp"], 0),
  ymax = c(pt["B_ee"], pt["B_ee"]+pt["sh_exp"], pt["B_ee"]+pt["sh_exp"]+pt["sh_pop"], pt["B_ll"]),
  fill = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"))
# 误差棒：Early/Late 用 B_ee/B_ll 区间；增量段用 Shapley 分量区间（叠加在段顶）
eb <- data.frame(
  x = 1:4,
  y = c(pt["B_ee"], pt["B_ee"]+pt["sh_exp"], pt["B_ee"]+pt["sh_exp"]+pt["sh_pop"], pt["B_ll"]),
  lo = c(quantile(mc[,"B_ee"],0.025), pt["B_ee"]+quantile(mc[,"sh_exp"],0.025),
         pt["B_ee"]+pt["sh_exp"]+quantile(mc[,"sh_pop"],0.025), quantile(mc[,"B_ll"],0.025)),
  hi = c(quantile(mc[,"B_ee"],0.975), pt["B_ee"]+quantile(mc[,"sh_exp"],0.975),
         pt["B_ee"]+pt["sh_exp"]+quantile(mc[,"sh_pop"],0.975), quantile(mc[,"B_ll"],0.975)))

pa <- ggplot(wf) +
  geom_rect(aes(xmin=x-0.36, xmax=x+0.36, ymin=ymin, ymax=ymax, fill=fill), color=NA) +
  geom_segment(data=data.frame(x=1:3, y=c(pt["B_ee"], pt["B_ee"]+pt["sh_exp"], pt["B_ll"])),
               aes(x=x+0.36, xend=x+1-0.36, y=y, yend=y), linetype="dashed", color="#718096", linewidth=0.4) +
  geom_errorbar(data=eb, aes(x=x, ymin=lo, ymax=hi), width=0.14, color="#2D3748", linewidth=0.5) +
  geom_text(aes(x=x, y=ymax, label=c("1,211","+3,663","+1,210","6,084")),
            vjust=-1.6, size=2.7, color="#2D3748", fontface="bold") +
  scale_fill_identity() +
  scale_x_continuous(breaks=1:4, labels=wf$label) +
  scale_y_continuous(labels=scales::comma, expand=expansion(mult=c(0,0.12))) +
  labs(tag="a", x=NULL, y="Modeled burden (thousand YLDs)") +
  theme_nature()

# ── 面板 b：水平贡献条 + 95% 区间 ───────────────────────────────────────────
share <- data.frame(
  y = c(2, 1),
  label = c("Warm-year exposure", "Female population ≥55"),
  est = c(pt["exp_share"], pt["pop_share"]) * 100,
  lo = c(quantile(mc[,"exp_share"],0.025), quantile(mc[,"pop_share"],0.025)) * 100,
  hi = c(quantile(mc[,"exp_share"],0.975), quantile(mc[,"pop_share"],0.975)) * 100,
  fill = c("#ff7f0e", "#2ca02c"))

pb <- ggplot(share, aes(y=reorder(label, -y), x=est)) +
  geom_col(aes(fill=fill), width=0.52) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.14, color="#2D3748", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.1f%%", est)), hjust=-0.15, size=3, color="#2D3748", fontface="bold") +
  scale_fill_identity() +
  scale_x_continuous(limits=c(0,100), breaks=seq(0,100,25), labels=function(x) paste0(x,"%"),
                     expand=expansion(mult=c(0,0.02))) +
  labs(tag="b", x="Modeled contribution share", y=NULL) +
  theme_nature() +
  theme(axis.text.y=element_text(color="#2D3748", size=8.5))

fig <- pa + pb + plot_layout(widths=c(1.45, 1))
ggsave("Figure4_main.png", fig, width=183, height=95, units="mm", dpi=300, bg="white")
ggsave("Figure4_main.pdf", fig, width=183, height=95, units="mm", device=cairo_pdf, bg="white")
cat("Fig 4 v9 saved\n")
