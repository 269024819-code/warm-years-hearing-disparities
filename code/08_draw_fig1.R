#!/usr/bin/env Rscript
# Fig 1 v9（老师意见）：a=点估计+误差线+β标注+显著性+0线；b=简化时间序列（注明固定5 Hot年）
library(ggplot2); library(dplyr); library(patchwork)

theme_nature <- function(base_size=9) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          panel.grid.major.y=element_line(color="#E2E8F0", linewidth=0.4),
          axis.line=element_line(color="#4A5568", linewidth=0.4),
          axis.ticks=element_line(color="#4A5568", linewidth=0.3),
          plot.title=element_text(face="bold", size=11, hjust=0),
          axis.title=element_text(color="#2D3748", size=9),
          axis.text=element_text(color="#4A5568", size=8),
          plot.tag=element_text(face="bold", size=12),
          text=element_text(family="Liberation Sans"))
}

# ── 面板a：Table 1 原始值（DK CI）──────────────────────────────────────────
t1 <- read.csv("data/table1_modsevere_160.csv")
t1$age_mid <- c(57,62,67,72,77,82)
t1$sig <- ifelse(t1$sd_p < 0.001, "***", ifelse(t1$sd_p < 0.01, "**", ifelse(t1$sd_p < 0.05, "*", "")))
t1$label <- sprintf("%.2f%s", t1$sexdiff_b, t1$sig)

pa <- ggplot(t1, aes(x=age_mid, y=sexdiff_b)) +
  geom_hline(yintercept=0, color="#4A5568", linewidth=0.6, linetype="solid") +
  geom_errorbar(aes(ymin=sd_lo, ymax=sd_hi), width=1.6, color="#1f77b4", linewidth=0.7) +
  geom_point(color="#1f77b4", size=2.8) +
  geom_text(aes(label=label, y=sd_hi), vjust=-0.7, size=2.8, color="#2D3748") +
  scale_x_continuous(breaks=t1$age_mid, labels=c("55-59","60-64","65-69","70-74","75-79","80+")) +
  labs(tag="a", x="Age group (years)",
       y="Hot-year × female interaction (YLDs per 100,000)") +
  coord_cartesian(ylim=c(-0.6, 7.4)) +
  theme_nature()

# ── 面板b：Hot年占比时间序列（简化）─────────────────────────────────────────
tp <- read.csv("data/hearing_temperature_analysis_panel_2003_2019.csv")
st <- readRDS("data/panel_stacked.rds")
temp <- tp %>% distinct(iso3, year, .keep_all=TRUE) %>%
  dplyr::select(iso3, year, temp_c=ts_annual_mean_temp_c) %>% filter(iso3 %in% unique(st$iso3))
hot <- temp %>% group_by(iso3) %>% mutate(hot=as.integer(temp_c>=quantile(temp_c,0.75,type=1))) %>% ungroup()
ts <- hot %>% group_by(year) %>% summarise(hot_share=mean(hot), .groups="drop")

pb <- ggplot(ts, aes(x=year, y=hot_share)) +
  geom_col(fill="#ff7f0e", width=0.72) +
  scale_x_continuous(breaks=seq(2003,2019,4)) +
  scale_y_continuous(labels=function(x) paste0(round(x*100),"%"), expand=expansion(mult=c(0,0.08))) +
  labs(tag="b", x="Year", y="Share of 160 countries in a hot year") +
  theme_nature()

fig <- pa / pb + plot_layout(heights=c(1.15, 1))
ggsave("Figure1_main.png", fig, width=183, height=150, units="mm", dpi=300, bg="white")
ggsave("Figure1_main.pdf", fig, width=183, height=150, units="mm", device=cairo_pdf, bg="white")
cat("Fig 1 v9 saved | 1b均值:", round(mean(ts$hot_share),3), "\n")
