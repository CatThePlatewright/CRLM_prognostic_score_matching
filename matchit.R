library(cobalt)
library(randomForest)
library("MatchIt")


train_data <- read.csv("traindata/unmatched_train_data.csv", stringsAsFactors = FALSE)
train_data$five_y_death <- ifelse(train_data$died == 1 & train_data$time_to_event <= 365*5, 1, 0)

# No matching; constructing a pre-match matchit object
#m.out0 <- matchit(adjuvantchemo~ age+T_stage+N_stage+
#                  rightleft.Rec + cea_carcinoembryonic.antigen+ extrahepatic.disease+
  #                size+ number_liver_mets+ R0+KRAS, data = train_data,
   #               method = NULL, distance = "glm")


tutorial<-train_data

tutorial$adjuvantchemo <- factor(tutorial$adjuvantchemo)

rf.out <- randomForest(adjuvantchemo ~ age+T_stage+N_stage+
                         rightleft.Rec + cea_carcinoembryonic.antigen+
                         extrahepatic.disease+
                         size+ number_liver_mets+ R0+KRAS,
                data=tutorial)

rf.out

eps <- rf.out$votes[,2] # Estimated PS 


match_obj <- matchit(adjuvantchemo ~ age+T_stage+N_stage+
                       rightleft.Rec + cea_carcinoembryonic.antigen+ extrahepatic.disease+
                       size+ number_liver_mets+ R0+KRAS, data=tutorial,
                     
                     distance=eps, # our own estimated PS
                     
                     method='nearest',
                     
                     replace=FALSE,
                     
                     discard='none',
                     
                     ratio=1)

match_obj <- matchit(adjuvantchemo ~ age+T_stage+N_stage+
                       rightleft.Rec + cea_carcinoembryonic.antigen+ extrahepatic.disease+
                       size+ number_liver_mets+ R0+KRAS, data = tutorial, method = "nearest",
                     distance = "glm", link = "linear.logit",
                     caliper = .2)

# Extract summary details as a dataframe
summary_obj <- summary(match_obj, standardize = TRUE, interactions = TRUE, un=FALSE)

# Save balance table (standardized mean differences, etc.)
balance_table <- summary_obj$sum.matched  # Extract balance statistics for matched samples
write.csv(balance_table, "matchit_glm_balance_table.csv", row.names = FALSE)

match_obj$match.matrix


plot(match_obj, type = "density", interactive = FALSE,
     
     which.xs = ~age+T_stage+N_stage+
       rightleft.Rec + cea_carcinoembryonic.antigen+ extrahepatic.disease+
       size+ number_liver_mets+ R0+KRAS)



m.data <- match.data(match_obj)

head(m.data)



write.csv(m.data,"traindata/matchit_train_glm.csv")

bal.tab(match_obj, m.threshold = 0.1, un = TRUE)
# Generate and save balance plots
plot <- bal.plot(match_obj, var.name = 'age', which = 'both', grid=TRUE) + ggtitle("age")
ggsave(filename = "balance_plot_age.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'T_stage', which = 'both', grid=TRUE) + ggtitle("T_stage")
ggsave(filename = "balance_plot_T_stage.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'N_stage', which = 'both', grid=TRUE) + ggtitle("N_stage")
ggsave(filename = "balance_plot_N_stage.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'rightleft.Rec', which = 'both', grid=TRUE) + ggtitle("rightleft.Rec")
ggsave(filename = "balance_plot_rightleft_Rec.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'cea_carcinoembryonic.antigen', which = 'both', grid=TRUE) + ggtitle("cea_carcinoembryonic.antigen")
ggsave(filename = "balance_plot_cea_carcinoembryonic_antigen.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'extrahepatic.disease', which = 'both', grid=TRUE) + ggtitle("extrahepatic.disease")
ggsave(filename = "balance_plot_extrahepatic_disease.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'size', which = 'both', grid=TRUE) + ggtitle("size")
ggsave(filename = "balance_plot_size.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'number_liver_mets', which = 'both', grid=TRUE) + ggtitle("number_liver_mets")
ggsave(filename = "balance_plot_number_liver_mets.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'R0', which = 'both', grid=TRUE) + ggtitle("R0")
ggsave(filename = "balance_plot_R0.png", plot = plot, width = 8, height = 6)

plot <- bal.plot(match_obj, var.name = 'KRAS', which = 'both', grid=TRUE) + ggtitle("KRAS")
ggsave(filename = "balance_plot_KRAS.png", plot = plot, width = 8, height = 6)

love.plot(bal.tab(match_obj, m.threshold=0.1),
          stat = "mean.diffs",
          grid=TRUE,
          stars="raw",
          abs = F)

