
library("survminer")
library("survival")
library("tidyverse")
library("Rcpp")


data <- read.csv("traindata/unmatched_train_data.csv")
#data <- read.csv("traindata/matched_train_equalized_r1.csv")
untreated_train_data_3a <- subset(data, adjuvantchemo == 0)
treated_train_data_3b <- subset(data, adjuvantchemo == 1)

titlename <-"KM plot of Cohort before Matching"
pngname <- "plots_and_metrics/KMplot_unmatched.png"
tutorial <- treated_train_data_3b
tutorial <- untreated_train_data_3a
# Convert time_to_event from days to months
tutorial$time_to_event <- tutorial$time_to_event / 30.44
library(survival)
library(prodlim)
quantile(prodlim(Surv(time=time_to_event, event=died)~1, data=tutorial, reverse=T))
library("gtsummary")
library("survminer")
survfit(Surv(time_to_event, died) ~ 1, data = tutorial) %>%
  tbl_survfit(
    times = c(60),
    label_header = "**5-year survival (95% CI)**"
  )
surv_obj <- survfit(Surv(time_to_event, died) ~ 1, data = tutorial)

# Extract median survival time
median_surv <- surv_median(surv_obj)

# Print the median survival time
print(paste("Median Overall Survival: ", median_surv$median, "days"))
print(paste("95% CI: ", median_surv$lower, " - ", median_surv$upper, "days"))


fit <- survfit(Surv(time_to_event, died)~adjuvantchemo, data = tutorial)
# Define a color palette with 10 distinct colors
risk_palette <- c("red", "green", "blue", "orange", "purple", "pink", "brown", "gray", "cyan", "magenta")

# Create Kaplan-Meier plot with different colors for risk buckets
KM_plot <- ggsurvplot(fit, data = tutorial,
                      xlim = c(0, 60),
                      pval = TRUE,
                      risk.table = TRUE,
                      tables.theme = theme_cleantable(),
                      tables.height = 0.3,
                      fontsize = 5,
                      legend.title = element_blank(),
                      break.time.by = 6,
                      xlab = "Overall survival, months",
                      censor = FALSE,
                      font.x = c(12, "bold"),
                      ylab = "Survival probability",
                      font.y = c(12, "bold"),
                      font.tickslab = c(14, "plain"),
                      conf.int = TRUE,
                      pval.size = 5,
                      break.y.by = 0.2,
                      pval.coord = c(0, 0.03),
                      palette = risk_palette,  
                      legend = c(0.7, 0.9),
                      title = titlename
                      )  # Use different linetypes for risk buckets

# Save the plot
ggsave(pngname, plot = KM_plot$plot)
print(KM_plot)
