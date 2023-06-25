
library(tidyverse)
library(lubridate)
library(scales)

setwd("~/OneDrive - Kemri Wellcome Trust/Visualization/Cemia_Dash/data/")
data = read.csv("coverage_dataset.csv")

data$Date = as.Date(data$Date)

png("../assets/coverage_plot.png", width = 600,height = 300)
ggplot(data = data, aes(x = Date,y=average_cov)) +
  geom_line(aes(y = average_cov),color="firebrick")+
  geom_point(aes(y=average_cov),color="firebrick")+
  geom_linerange(aes(ymin = min_cov,ymax=max_cov),color="blue")+
  scale_x_date( breaks = date_breaks("3 months"),minor_breaks = "1 month",  labels = date_format("%b-%Y"))+
  scale_y_continuous(limits = c(0,32000), breaks = seq(0,32000,by=6000), labels = unit_format(unit = "K",scale=1e-3))+
  theme(axis.text.x = element_text(hjust=0.5,size=8,color="black"),axis.line.y = element_line(color="black", size=0.3),
        axis.title.x = element_blank(),panel.grid.major.y = element_line(colour = "gray75",size=0.1),
        panel.grid.major.x = element_line(colour = "gray75",size=0.1),axis.line.x = element_line(color="black", size=0.3),
        axis.text.y = element_text(size=8,color="black"),panel.background = element_blank())+
  labs(y = "Average Monthly Coverage")

dev.off()
