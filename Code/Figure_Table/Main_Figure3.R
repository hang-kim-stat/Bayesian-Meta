# Load necessary libraries
library(ggplot2)
library(patchwork) # For arranging plots side by side
library(dplyr)
library(tidyr)
library(scales)

# Set text size
SIZE_title.y = 18 ; SIZE_title.x = 16 ; SIZE_text.y = 16 ; SIZE_text.x = 18 ; SIZE_title = 22 ; SIZE_b = 10

# Define Greek labels for the legend
greek_labels <- c(expression(mu[0]), expression(mu[1]), expression(mu[2]), expression(mu[12]))

##############################################################################
# Load the results of 1st simulation 

# data1 = read.csv("../../output/Figure_Table/Supplement_Table4.csv", stringsAsFactors = FALSE)
Table4 = read.csv("../../output/Figure_Table/Supplement_Table4.csv")

Table4_figure = array(NA,c(16,5))
colnames(Table4_figure) = c("Method","Parameter","Bias","MSE","Coverage")
Table4_figure[,1] = rep(c("Benchmark","IPD-AD","IPD-AD(pooled)","IPD only"),4)
Table4_figure[,2] = rep(c("mu1","mu2","mu3","mu4"),each=4)
for (i_method in 1:4){
 SEQ1 = c(0,4,8,12)+i_method
 SEQ2 = (2:4)+(i_method-1)*4
 Table4_figure[SEQ1,3:5] = as.matrix(Table4[1:4,SEQ2])
} # 

# Set the order of the Method variable
Table4_figure[,"Method"] <- factor(Table4_figure[,"Method"], levels = c("Benchmark", "IPD-AD", "IPD-AD(pooled)", "IPD only"), labels = c("Bench\nmark", "IPD-AD", "IPD-AD\n(pooled)", "IPD\nonly") )

# Plot for Bias with legend
plot_bias1 <- ggplot(Table4_figure, aes(x = Method, y = Bias, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) +
 geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
 labs(title = "Bias", x = "", y = "Bias") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(-0.1, 0.1), labels = number_format(accuracy = 0.01)) +  
 theme(
  plot.title = element_text(size = SIZE_title, face = "bold", hjust = 0.5, margin = margin(b = SIZE_b)),
  axis.title.x = element_blank(),
  # axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 20)),
  # axis.title.y = element_text(size = SIZE_title.y),
  axis.title.y = element_blank(),
  # axis.text.x = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.position = "none"
 )

# Plot for MSE with connecting lines, without legend, and with adjusted text size
plot_mse1 <- ggplot(Table4_figure, aes(x = Method, y = MSE, group = Parameter, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) + 
 geom_line(linetype = "dashed") +
 labs(title = "MSE", x = "", y = "MSE") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(0, 0.13)) +  
 theme(
  plot.title = element_text(size = SIZE_title, face = "bold", hjust = 0.5, margin = margin(b = SIZE_b)),
  legend.position = "none",
  axis.title.x = element_blank(),
  # axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 20)),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.text = element_text(size = 16)    
 )

# Plot for 95% Coverage without legend and with adjusted text size
plot_95cover1 <- ggplot(Table4_figure, aes(x = Method, y = Coverage, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) +
 geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
 labs(title = "95% C.I. Coverage", x = "", y = "95% C.I. Coverage") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(0.9, 1.0), labels = number_format(accuracy = 0.01)) +  
 theme(
  plot.title = element_text(size = SIZE_title, face = "bold", hjust = 0.5, margin = margin(b = SIZE_b)),
  axis.title.x = element_blank(),
  # axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 20)),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.position = "none" 
 )

##############################################################################
# Load the data from CSV of 2nd simulation 

Table5 = read.csv("../../output/Figure_Table/Supplement_Table5.csv")

data2 <- read.csv("W_Simul2_Result_for_PaperIn_ChatGPT1.csv", stringsAsFactors = FALSE)

# Set the order of the Method variable
data2$Method <- factor(data2$Method, levels = c("Benchmark", "IPD-AD", "IPD-AD(pooled)", "IPD only"), labels = c("Bench\nmark", "IPD-AD", "IPD-AD\n(pooled)", "IPD\nonly") )

# Plot for Bias 
plot_bias2 <- ggplot(data2, aes(x = Method, y = Bias, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) +
 geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
 labs(x = "", y = "Bias") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(-0.1, 0.1), labels = number_format(accuracy = 0.01)) + 
 theme(
  # plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
  legend.position = "none", 
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.text = element_text(size = 16)    
 )

# Plot for MSE with connecting lines, without legend, and with adjusted text size
plot_mse2 <- ggplot(data2, aes(x = Method, y = MSE, group = Parameter, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) + 
 geom_line(linetype = "dashed") +
 labs(title = "", x = "", y = "MSE") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(0, 0.018), labels = number_format(accuracy = 0.01)) +  
 theme(
  # plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
  legend.position = "none",
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.text = element_text(size = 16)   
 )

# Plot for 95% Coverage without legend and with adjusted text size
plot_95cover2 <- ggplot(data2, aes(x = Method, y = Coverage, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) +
 geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
 labs(x = "", y = "95% C.I. Coverage") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(0.87, 1.0)) + 
 theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.position = "none"
 )

# # Combine the plots side by side
# combined_plot <- plot_bias1 + plot_bias2 + plot_mse1 + plot_mse2 + plot_95cover1  + plot_95cover2 +  plot_layout(nrow = 3, ncol = 2)
# 
# # Save the combined plot as a PNG file
# ggsave("Simul12_combined_v2024_0109.png", plot = combined_plot, width = 16, height = 16, dpi = 300)

# Combine the plots side by side
combined_plot <- plot_bias1 + plot_mse1 + plot_95cover1 + plot_bias2 + plot_mse2 +  plot_95cover2 + plot_layout(nrow = 2, ncol = 3)

ggsave("Simul12_combined_v2024_0208.png", plot = combined_plot, width = 16, height = 8, dpi = 300)

# Plot for 95% Coverage without legend and with adjusted text size
plot_95cover2_for_legend <- ggplot(data2, aes(x = Method, y = Coverage, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) +
 geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
 labs(x = "", y = "") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(0.85, 1.0)) + 
 theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.position = "bottom"
 )

ggsave("For_legends.png", plot = plot_95cover2_for_legend, width = 8, height = 6, dpi = 300)
