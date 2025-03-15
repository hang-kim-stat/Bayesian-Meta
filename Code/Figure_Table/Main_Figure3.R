rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(patchwork) # For arranging plots side by side
library(dplyr)
library(tidyr)
library(scales)

# Set text size
SIZE_title.y = 18 ; SIZE_title.x = 16 ; SIZE_text.y = 16 ; SIZE_text.x = 18 ; SIZE_title = 22 ; SIZE_b = 10

# Define greek labels for the legend
greek_labels <- c(expression(mu[0]), expression(mu[1]), expression(mu[2]), expression(mu[12]))

##############################################################################

# Load the results of 1st simulation 
Table4 = read.csv("../../Output/Figure_Table/Supplement_Table4.csv")

Table4_label = array(NA,c(16,2))
colnames(Table4_label) = c("Method","Parameter")
Table4_label[,1] = rep(c("Benchmark","IPD-AD","IPD-AD(pooled)","IPD only"),4)
Table4_label[,2] = rep(c("mu1","mu2","mu3","mu4"),each=4)

Table4_number = array(NA,c(16,3))
colnames(Table4_number) = c("Bias","MSE","Coverage")
for (i_method in 1:4){
 SEQ1 = c(0,4,8,12)+i_method
 SEQ2 = (3:5)+(i_method-1)*4
 tempDataFrame = Table4[1:4,SEQ2]
 tempVector = as.numeric(as.matrix(tempDataFrame))
 Table4_number[SEQ1,] = matrix(tempVector, nrow = nrow(tempDataFrame), ncol = ncol(tempDataFrame))
}

Table4_figure <- data.frame(Method = Table4_label[,"Method"], Parameter = Table4_label[,"Parameter"], Bias = Table4_number[,"Bias"], MSE = Table4_number[,"MSE"], Coverage = Table4_number[,"Coverage"], stringsAsFactors = FALSE)

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
  axis.title.y = element_blank(),
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

# Load the results of 2nd simulation 
Table5 = read.csv("../../Output/Figure_Table/Supplement_Table5.csv")

Table5_label = array(NA,c(16,2))
colnames(Table5_label) = c("Method","Parameter")
Table5_label[,1] = rep(c("Benchmark","IPD-AD","IPD-AD(pooled)","IPD only"),4)
Table5_label[,2] = rep(c("mu1","mu2","mu3","mu4"),each=4)

Table5_number = array(NA,c(16,3))
colnames(Table5_number) = c("Bias","MSE","Coverage")
for (i_method in 1:4){
 SEQ1 = c(0,4,8,12)+i_method
 SEQ2 = (3:5)+(i_method-1)*4
 tempDataFrame = Table5[1:4,SEQ2]
 tempVector = as.numeric(as.matrix(tempDataFrame))
 Table5_number[SEQ1,] = matrix(tempVector, nrow = nrow(tempDataFrame), ncol = ncol(tempDataFrame))
}

Table5_figure <- data.frame(Method = Table5_label[,"Method"], Parameter = Table5_label[,"Parameter"], Bias = Table5_number[,"Bias"], MSE = Table5_number[,"MSE"], Coverage = Table5_number[,"Coverage"], stringsAsFactors = FALSE)

# Set the order of the Method variable
Table5_figure[,"Method"] <- factor(Table5_figure[,"Method"], levels = c("Benchmark", "IPD-AD", "IPD-AD(pooled)", "IPD only"), labels = c("Bench\nmark", "IPD-AD", "IPD-AD\n(pooled)", "IPD\nonly") )

# Plot for Bias 
plot_bias2 <- ggplot(Table5_figure, aes(x = Method, y = Bias, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) +
 geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
 labs(x = "", y = "Bias") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(-0.1, 0.1), labels = number_format(accuracy = 0.01)) + 
 theme(
  legend.position = "none", 
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.text = element_text(size = 16)    
 )

# Plot for MSE with connecting lines, without legend, and with adjusted text size
plot_mse2 <- ggplot(Table5_figure, aes(x = Method, y = MSE, group = Parameter, shape = Parameter, color = Parameter)) +
 geom_point(size = 5) + 
 geom_line(linetype = "dashed") +
 labs(title = "", x = "", y = "MSE") +
 scale_shape_manual(values = c(1, 2, 15, 4), labels = greek_labels) +
 scale_color_manual(values = c("blue", "brown", "red", "purple"), labels = greek_labels) +
 scale_y_continuous(limits = c(0, 0.018), labels = number_format(accuracy = 0.01)) +  
 theme(
  legend.position = "none",
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = SIZE_text.x),
  axis.text.y = element_text(size = SIZE_text.y),
  legend.text = element_text(size = 16)   
 )

# Plot for 95% Coverage without legend and with adjusted text size
plot_95cover2 <- ggplot(Table5_figure, aes(x = Method, y = Coverage, shape = Parameter, color = Parameter)) +
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

##############################################################################

# Combine the plots side by side
combined_plot <- plot_bias1 + plot_mse1 + plot_95cover1 + plot_bias2 + plot_mse2 +  plot_95cover2 + plot_layout(nrow = 2, ncol = 3)

# Store figure in .png file
ggsave("../../Output/Figure_Table/Maintext_Figure3.png", plot = combined_plot, width = 16, height = 8, dpi = 300)