#SCOP flow cytometry data
####Set up environment####
getwd()
setwd("C:/SCOP/")
setwd("G:/SCOP/")

library(plyr)
library(tidyverse)
library(ggpubr)

flow.cyto <- read.csv("SCOP_FlowCytoData.csv", header = T, sep = ",")

####Environmental data####
SCOP_meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",")
SCOP_meta$Date_Combined <- paste(SCOP_meta$Correct_Month, SCOP_meta$Correct_Day, SCOP_meta$Correct_Year, sep = "/")
SCOP_meta$Correct_Date <- as.Date(SCOP_meta$Date_Combined, "%m/%d/%Y")

####Microscopy data####
SCOP.phyto <- read.csv("SCOP_PhytoplanktonCounts_Avg2.csv", header = T, sep = ",")
SCOP.phyto <- SCOP.phyto[, -c(1)]
SCOP.phyto$Date <- as.Date(SCOP.phyto$Date_Long, format = "%d-%b-%y")
SCOP.phyto2 <- SCOP.phyto[, c(3, 4, 7)]
SCOP.phyto.spread <- SCOP.phyto2 %>% spread(PlanktonGroup, Average_Cells_mL)

####Synechococcus####
flow.cyto.Syn <- ddply(flow.cyto,~Date,summarise,mean=mean(Syn.Cell.Density..cells.mL.),sd=sd(Syn.Cell.Density..cells.mL.))
colnames(flow.cyto.Syn) <- c("Date", "Mean.Syn.CellDensity", "SD.Syn.CellDensity")
flow.cyto.Syn$Date <- c("10/11/2021", "10/13/2021", "10/15/2021", "10/18/2021", "10/20/2021", "10/22/2021", "10/25/2021", "10/27/2021", "10/29/2021", "10/06/2021", "10/08/2021", "11/10/2021", "11/17/2021", "11/24/2021", "11/03/2021", "12/01/2021", "12/15/2021", "12/22/2021", "12/29/2021", "12/08/2021")
flow.cyto.Syn$Date <- as.character(flow.cyto.Syn$Date)
flow.cyto.Syn$Date <- as.Date(flow.cyto.Syn$Date, "%m/%d/%Y")
flow.cyto.Syn$Timeline <- c(rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Week 1", 2), rep("Month 2", 4), rep("Month 3", 5))

range(flow.cyto.Syn$Mean.Syn.CellDensity) #21231.10 95686.65
syn.plot <- ggplot() + 
  geom_line(data = flow.cyto.Syn, aes(x = Date, y = Mean.Syn.CellDensity), linetype = "dotted", linewidth = 0.5, color = "black") + 
  geom_point(data = flow.cyto.Syn, aes(x = Date, y = Mean.Syn.CellDensity, color = Timeline), shape = 16, size = 3) + 
  theme_classic() + 
  geom_errorbar(data = flow.cyto.Syn, aes(x = Date, y = Mean.Syn.CellDensity, ymin = Mean.Syn.CellDensity-SD.Syn.CellDensity, ymax = Mean.Syn.CellDensity+SD.Syn.CellDensity), width= 1.5, position = position_dodge(0.05)) + 
  labs(x = "Date", y = paste("Avg. Synechococcus\n Cell Density (cells/mL)")) + 
  expand_limits(y = c(20000, 100000)) + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_date(breaks = flow.cyto.Syn$Date, date_labels = "%Y-%m-%d") + 
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c("transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", "grey30", "transparent", "grey30", "grey30", "grey30", "grey30", "grey30", "navy", "grey30", "firebrick3", "grey30"))) + 
  scale_color_manual(values = c("firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30"), name = "Timeline:", breaks = c("Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3"), labels = c("October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December"))
syn.plot

####Eukaryotes####
flow.cyto.Euks <- ddply(flow.cyto,~Date,summarise,mean=mean(Total.Euk.Cell.Density..cells.mL.),sd=sd(Total.Euk.Cell.Density..cells.mL.))
colnames(flow.cyto.Euks) <- c("Date", "Mean.Euks.CellDensity", "SD.Euks.CellDensity")
flow.cyto.Euks$Date <- c("10/11/2021", "10/13/2021", "10/15/2021", "10/18/2021", "10/20/2021", "10/22/2021", "10/25/2021", "10/27/2021", "10/29/2021", "10/06/2021", "10/08/2021", "11/10/2021", "11/17/2021", "11/24/2021", "11/03/2021", "12/01/2021", "12/15/2021", "12/22/2021", "12/29/2021", "12/08/2021")
flow.cyto.Euks$Date <- as.character(flow.cyto.Euks$Date)
flow.cyto.Euks$Date <- as.Date(flow.cyto.Euks$Date, "%m/%d/%Y")
flow.cyto.Euks$Timeline <- c(rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Week 1", 2), rep("Month 2", 4), rep("Month 3", 5))

range(flow.cyto.Euks$Mean.Euks.CellDensity) #4430.00 20206.67
euk.plot <- ggplot() + 
  geom_line(data = flow.cyto.Euks, aes(x = Date, y = Mean.Euks.CellDensity), linetype = "dotted", linewidth = 0.5, color = "black") + 
  geom_point(data = flow.cyto.Euks, aes(x = Date, y = Mean.Euks.CellDensity, color = Timeline), shape = 16, size = 3) + 
  geom_errorbar(data = flow.cyto.Euks, aes(x = Date, y = Mean.Euks.CellDensity, ymin = Mean.Euks.CellDensity-SD.Euks.CellDensity, ymax = Mean.Euks.CellDensity+SD.Euks.CellDensity), width= 1.5, position = position_dodge(0.05)) + 
  theme_classic() + 
  labs(x = "Date", y = "Avg. Eukaryote\n Cell Density (cells/mL)") + 
  expand_limits(y = c(5000, 21000)) + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_date(breaks = flow.cyto.Euks$Date, date_labels = "%Y-%m-%d") + 
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour =c("transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", "grey30", "transparent", "grey30", "grey30", "grey30", "grey30", "grey30", "navy", "grey30", "firebrick3", "grey30"))) + 
  scale_color_manual(values = c("firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30"), name = "Timeline:", breaks = c("Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3"), labels = c("October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December"))
euk.plot

####Heterotrophs####
flow.cyto.Hetero <- ddply(flow.cyto,~Date,summarise,mean=mean(Hetero.Cell.Density..cells.mL.),sd=sd(Hetero.Cell.Density..cells.mL.))
colnames(flow.cyto.Hetero) <- c("Date", "Mean.Hetero.CellDensity", "SD.Hetero.CellDensity")
flow.cyto.Hetero$Date <- c("10/11/2021", "10/13/2021", "10/15/2021", "10/18/2021", "10/20/2021", "10/22/2021", "10/25/2021", "10/27/2021", "10/29/2021", "10/06/2021", "10/08/2021", "11/10/2021", "11/17/2021", "11/24/2021", "11/03/2021", "12/01/2021", "12/15/2021", "12/22/2021", "12/29/2021", "12/08/2021")
flow.cyto.Hetero$Date <- as.character(flow.cyto.Hetero$Date)
flow.cyto.Hetero$Date <- as.Date(flow.cyto.Hetero$Date, "%m/%d/%Y")
flow.cyto.Hetero$Timeline <- c(rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Week 1", 2), rep("Month 2", 4), rep("Month 3", 5))

range(flow.cyto.Hetero$Mean.Hetero.CellDensity) #161330.0 949461.7
hetero.plot <- ggplot() + 
  geom_point(data = flow.cyto.Hetero, aes(x = Date, y = Mean.Hetero.CellDensity, color = Timeline), shape = 16, size = 3) + 
  geom_line(data = flow.cyto.Hetero, aes(x = Date, y = Mean.Hetero.CellDensity), linetype = "dotted", linewidth = 0.5, color = "black") + 
  geom_errorbar(data = flow.cyto.Hetero, aes(x = Date, y = Mean.Hetero.CellDensity, ymin = Mean.Hetero.CellDensity-SD.Hetero.CellDensity, ymax = Mean.Hetero.CellDensity+SD.Hetero.CellDensity), width= 1.5, position = position_dodge(0.05)) + 
  theme_classic() + 
  labs(x = "Date", y = "Avg. Heterotrophic Bacterial\n Cell Density (cells/mL)") + 
  expand_limits(y = c(160000, 1000000)) + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_date(breaks = flow.cyto.Hetero$Date, date_labels = "%Y-%m-%d") + 
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c("transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", "grey30", "transparent", "grey30", "grey30", "grey30", "grey30", "grey30", "navy", "grey30", "firebrick3", "grey30"))) + 
  scale_color_manual(values = c("firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30"), name = "Timeline:", breaks = c("Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3"), labels = c("October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December"))
hetero.plot

####Phytoplankton - Microscopy####
#With flagellates divided by 100
SCOP.phyto100 <- read.csv("SCOP_PhytoplanktonCounts_Avg2.csv", header = T, sep = ",")
SCOP.phyto100$Date_Long <- as.Date(SCOP.phyto100$Date_Long, format = "%d-%b-%y")

SCOP.phyto.plot3 <- ggplot(data = SCOP.phyto100, mapping = aes(x = Date_Long, y = Average_Cells_mL, fill = factor(PlanktonGroup, levels = c("Centric Diatoms", "Pennate Diatoms", "Dinoflagellates", "Flagellates/100"))))+
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Average_Cells_mL-Standard_Dev, ymax = Average_Cells_mL+Standard_Dev), width = 0.6, colour="black", linewidth = 0.5, position = position_dodge(1.75)) + 
  theme_classic() + 
  labs(y = "Avg. Cell Density (cells/mL)", x = "Date", fill = "Plankton Group:") + 
  expand_limits(y = c(0, 450)) + 
  scale_y_continuous(breaks = seq(0, 450, 50)) + 
  scale_x_date(breaks = SCOP.phyto100$Date_Long, date_labels = "%Y-%m-%d") +
  scale_fill_manual(breaks = c("Centric Diatoms", "Pennate Diatoms", "Dinoflagellates", "Flagellates/100"), labels = c("Centric Diatoms", "Pennate Diatoms", "Dinoflagellates", "Flagellates/100"), values = c("darkseagreen3", "darkseagreen4", "lightsteelblue3", "plum4")) +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, color = c("grey30", "grey30", "grey30", "grey30", "seagreen4", "firebrick3", "grey30", "grey30", "grey30", "grey30", "grey30")))
SCOP.phyto.plot3

####Combined plots####
setwd("G:/SCOP/Figures/")
getwd()

fig4.p1 <- ggarrange(syn.plot, euk.plot, hetero.plot, ncol = 3, nrow = 1, align = "hv", common.legend = T, legend = "right", labels = c("b", "c", "d"), font.label = list(size = 14, color = "black", face = "bold", family = "sans"))

png("fig4v2.png", width = 12, height = 7, units = 'in', res = 600)
ggarrange(SCOP.phyto.plot3, fig4.p1, ncol = 1, nrow = 2, labels = c("a"), font.label = list(size = 14, color = "black", face = "bold", family = "sans"))
dev.off()
