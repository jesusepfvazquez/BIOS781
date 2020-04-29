# Code for facet wrapped box plots
# Get data
scenario1 = read.csv("m = 10_n = 500_p = 10000_corr.csv")
scenario2 = read.csv("m = 100_n = 500_p = 10000_corr.csv")
scenario3 = read.csv("m = 1000_n = 500_p = 10000_corr.csv")
scenario4 = read.csv("m = 10000_n = 500_p = 10000_corr.csv")
scenario5 = read.csv("m = 10_n = 1000_p = 10000_corr.csv")
scenario6 = read.csv("m = 100_n = 1000_p = 10000_corr.csv")
scenario7 = read.csv("m = 1000_n = 1000_p = 10000_corr.csv")
scenario8 = read.csv("m = 10000_n = 1000_p = 10000_corr.csv")
scenario9 = read.csv("m = 10_n = 5000_p = 10000_corr.csv")
scenario10 = read.csv("m = 100_n = 5000_p = 10000_corr.csv")
scenario11 = read.csv("m = 1000_n = 5000_p = 10000_corr.csv")
scenario12 = read.csv("m = 10000_n = 5000_p = 10000_corr.csv")


# Label data for later groupings
scenario1$label=("n/p=0.05, m/p=0.001")
scenario2$label=("n/p=0.05, m/p=0.01")
scenario3$label=("n/p=0.05, m/p=0.1")
scenario4$label=("n/p=0.05, m/p=1")
scenario5$label=("n/p=0.1, m/p=0.001")
scenario6$label=("n/p=0.1, m/p=0.01")
scenario7$label=("n/p=0.1, m/p=0.1")
scenario8$label=("n/p=0.1, m/p=1")
scenario9$label=("n/p=0.5, m/p=0.001")
scenario10$label=("n/p=0.5, m/p=0.01")
scenario11$label=("n/p=0.5, m/p=0.1")
scenario12$label=("n/p=0.5, m/p=1")

# stack all 12 datasets
alldat = rbind(scenario1, scenario2, scenario3, scenario4, scenario5, scenario6, scenario7, scenario8, scenario9, scenario10, scenario11, scenario12)

# Order the labels according to how we would like them presented in the final plot
alldat$label <- factor(alldat$label, levels = c("n/p=0.05, m/p=0.001", 
                                                "n/p=0.05, m/p=0.01", 
                                                "n/p=0.05, m/p=0.1", 
                                                "n/p=0.05, m/p=1",
                                                "n/p=0.1, m/p=0.001",
                                                "n/p=0.1, m/p=0.01", 
                                                "n/p=0.1, m/p=0.1", 
                                                "n/p=0.1, m/p=1",
                                                "n/p=0.5, m/p=0.001", 
                                                "n/p=0.5, m/p=0.01", 
                                                "n/p=0.5, m/p=0.1",
                                                "n/p=0.5, m/p=1"))


# Transpose data so p_cuts become column names for all X1:X100 values
alldat_t <- alldat %>%
  gather(key, value, X1:X100)

# Probably don't need to do this, just deleting extra column from the transpose
alldat_t$key = NULL

# Similar to above, make p_cuts a factor and apply p-values from largest to smallest to appear in plots
alldat_t$p_cuts <- factor(alldat_t$p_cuts, levels = c("1", "0.1", "0.01", "0.001", "1e-05", "1e-08"))

# Graphics
p <- ggplot(alldat_t, aes(factor(p_cuts), value, fill=factor(p_cuts)))
pg <- p + geom_boxplot() + facet_wrap(~label, ncol=4) + ylab("Correlation") +xlab("P-value Cutoff") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Prediction Accuracy by Scenario") +
  theme(legend.position = "none", strip.text = element_text(size = 6), plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red") 
