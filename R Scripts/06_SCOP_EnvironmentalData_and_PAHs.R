####SCOP Summary
####Set up environment####
#Set working directory
getwd()
setwd("G:/SCOP/")
getwd()

#Load libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(ggtext)

#Read in metadata file
SCOP_meta <- read.csv("NP SCOP Metadata.csv", header = T, sep = ",")

#Convert date to date format
SCOP_meta$Date_Combined <- paste(SCOP_meta$Correct_Month, SCOP_meta$Correct_Day, SCOP_meta$Correct_Year, sep = "/")
SCOP_meta$Correct_Date <- as.Date(SCOP_meta$Date_Combined, "%m/%d/%Y")

#Add in timeline variable
SCOP_meta$Timeline <- c(rep("Before", 6), rep("Week 1", 3), rep("Week 2", 3), rep("Week 3", 3), rep("Week 4", 3), rep("Month 2", 4), rep("Month 3", 5), rep("Month 4", 4))

####Plots of Environmental Data (Temperature, Salinity, CUTI)####
#Temperature
range(SCOP_meta$Temperature_.C., na.rm = T) #14.5399 19.8672
SCOP.temp.plot <- ggplot(data = SCOP_meta, aes(x = Correct_Date, y = Temperature_.C.)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = "Temperature (°C)", x = "Date") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(14, 20)) + 
  scale_y_continuous(breaks = seq(14, 20, 1)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.temp.plot


#Salinity
range(SCOP_meta$Salinity, na.rm = T) #30.2749 33.5486
SCOP.sal.plot <-  ggplot(data = SCOP_meta, aes(x = Correct_Date, y = Salinity)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = "Salinity", x = "Date") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(30, 34)) + 
  scale_y_continuous(breaks = seq(30, 34, 1)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.sal.plot


#Upwelling Index
range(SCOP_meta$CUTI_33N, na.rm = T) #-0.170  0.883
SCOP.CUTI33N.plot <- ggplot(data = SCOP_meta, aes(x = Correct_Date, y = CUTI_33N)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = "CUTI (33N)", x = "Date") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(-0.25, 1.6)) + 
  scale_y_continuous(breaks = seq(-0.25, 1.5, 0.25)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.CUTI33N.plot

####Nutrient Concentrations (Nitrate, SRP)####
#Nitrate Concentrations
range(SCOP_meta$NO3_.uM._Avg, na.rm = T) #0  1.755123
SCOP.nitrate.plot <-  ggplot(data = SCOP_meta, aes(x = Correct_Date, y = NO3_.uM._Avg)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = expression(paste("Avg. Nitrate (", mu, "M)")), x = "Date") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(0, 1.76)) + 
  scale_y_continuous(breaks = seq(0, 1.75, 0.25)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.nitrate.plot


#SRP Concentrations
range(SCOP_meta$SRP_.uM._Avg, na.rm = T) #0.06841718 0.69716444
SCOP.SRP.plot <- ggplot(data = SCOP_meta, aes(x = Correct_Date, y = SRP_.uM._Avg)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = expression(paste("Avg. SRP (", mu, "M)")), x = "Date") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(0, 0.80)) + 
  scale_y_continuous(breaks = seq(0, 0.80, 0.10)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.SRP.plot

####Wave direction####
#Read in files
buoy.data.2021 <- read.delim("46253h2021.txt", header = TRUE, sep = "")
colnames(buoy.data.2021)

buoy.data.2022 <- read.delim("46253h2022.txt", header = T, sep = "")
colnames(buoy.data.2022)

#Subset to mean wave direction
wave.data.2021 <- buoy.data.2021[-c(1), c(1:3, 12)]
wave.data.2022 <- buoy.data.2022[-c(1), c(1:3, 12)]

#Convert to date format
wave.data.2021$Date_Combined <- paste(wave.data.2021$MM, wave.data.2021$DD, wave.data.2021$X.YY, sep = "/")
wave.data.2021$Date <- as.Date(wave.data.2021$Date_Combined, "%m/%d/%Y")
wave.data.2021$MM <- as.numeric(wave.data.2021$MM)

wave.data.2022$Date_Combined <- paste(wave.data.2022$MM, wave.data.2022$DD, wave.data.2022$X.YY, sep = "/")
wave.data.2022$Date <- as.Date(wave.data.2022$Date_Combined, "%m/%d/%Y")
wave.data.2022$MM <- as.numeric(wave.data.2022$MM)

#Subset to August - December 
wave.data.2021.sub <- wave.data.2021[which(wave.data.2021$MM >= 8), c(4, 6)]

wave.data.2022.sub <- wave.data.2022[which(wave.data.2022$MM == 1), c(4, 6)]
View(wave.data.2022.sub)

#Calculate mean daily wave direction
wave.data.2021.sub$MWD <- as.numeric(wave.data.2021.sub$MWD)
wave.data.2021.sub2 <- aggregate(MWD ~ Date, wave.data.2021.sub, mean)

wave.data.2022.sub$MWD <- as.numeric(wave.data.2022.sub$MWD)
wave.data.2022.sub2 <- aggregate(MWD ~ Date, wave.data.2022.sub, mean)

#Combine dataframes
wave.data <- rbind(wave.data.2021.sub2, wave.data.2022.sub2)

#Combine with other metadata
wave.data$Correct_Date <- wave.data$Date
wave.data.env <- merge(SCOP_meta, wave.data, by = "Correct_Date")

#Plot
SCOP.mwd.plot <- ggplot() + 
  geom_line(data = wave.data, aes(x = Date, y = MWD), color = "black", linewidth = 0.5, linetype = "dotted") + 
  geom_point(data = wave.data, aes(x = Date, y = MWD), shape = 21, color = "black", fill = "white", size = 3) + 
  geom_point(data = wave.data.env, aes(x = Correct_Date, y = MWD, color = Timeline), size = 5) + 
  theme_classic() + 
  labs(y = "Daily Avg. Wave Direction") + 
  expand_limits(y = c(150, 300)) + 
  scale_y_continuous(breaks = seq(150, 300, 30)) + 
  scale_x_date(breaks = wave.data.env$Correct_Date, date_labels =  "%Y-%m-%d") +
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January")) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4))))
SCOP.mwd.plot

####Particulate Organic Carbon (POC)####
range(SCOP_meta$POC_.uM._Avg, na.rm = T) #10.95728 265.65595
SCOP.POC.plot <- ggplot(data = SCOP_meta, aes(x = Correct_Date, y = POC_.uM._Avg)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 3) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = expression(paste("Avg. POC (", mu, "M)")), x = "Date") + 
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(0, 275)) + 
  scale_y_continuous(breaks = seq(0, 275, 25)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.POC.plot

####Chlorophyll####
range(SCOP_meta$Chlorophyll_.ug.L., na.rm = T) #0.57826 39.96142
SCOP.chl.plot1 <- ggplot(data = SCOP_meta, aes(x = Correct_Date, y = Chlorophyll_.ug.L.)) + 
  geom_point(aes(color = Timeline), shape = 16, size = 5) + 
  geom_line(color = "black", linetype = "dotted", linewidth = 0.5) + 
  theme_classic() + 
  labs(y = expression(paste("Chlorophyll (", mu, "g/L)")), x = "Date") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 22), legend.title = element_text(size = 24), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c(rep("grey30", 6), "transparent", "grey30", "transparent", "transparent", "seagreen4", "transparent", "transparent", "firebrick3", "transparent", "transparent", "grey30", "transparent", rep("grey30", 6), "navy", "grey30", "firebrick3", rep("grey30", 4)))) +
  expand_limits(y = c(0, 40)) + 
  scale_y_continuous(breaks = seq(0, 40, 5)) + 
  scale_x_date(breaks = SCOP_meta$Correct_Date, date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c("darkseagreen4", "firebrick4", "tomato3", "salmon2", "peachpuff2", "grey80", "grey30", "black"), name = "Timeline:", breaks = c("Before", "Week 1", "Week 2", "Week 3", "Week 4", "Month 2", "Month 3", "Month 4"), labels = c("Prior to Oil Spill", "October - Week 1", "October - Week 2", "October - Week 3", "October - Week 4", "November", "December", "January"))
SCOP.chl.plot1


####Polycyclic Aromatic Hydrocarbons (PAHs)####
np.pah2 <- read.csv("OilSpill-quant-20220907-PAHs_withdups.csv", sep = ",")
np.pah2$Date <- as.Date(np.pah2$Date, format = "%m/%d/%Y")
np.pah2.grouped <- melt(np.pah2, id = "Date")
np.pah2.grouped.avg <- np.pah2.grouped %>% group_by(Date, variable) %>% summarise(Average=mean(value, na.rm = T)) 
np.pah2.grouped.sd <- np.pah2.grouped %>% group_by(Date, variable) %>% summarise(SD = sd(value, na.rm = T)) 

np.pah2.avg.spread <- as.data.frame(spread(np.pah2.grouped.avg, variable, Average))
np.pah2.avg.spread2 <- np.pah2.avg.spread[, -1]
row.names(np.pah2.avg.spread2) <- np.pah2.avg.spread$Date
np.pah2.avg.spread2[is.na(np.pah2.avg.spread2)] <- 0
np.pah2.avg.spread2 <- mutate_all(np.pah2.avg.spread2, function(x) as.numeric(as.character(x)))
np.pah2.avg.spread2$Total_PAH <- rowSums(np.pah2.avg.spread2, na.rm = T)
np.pah2.avg.spread2$Date <- row.names(np.pah2.avg.spread2)
np.pah2.avg.grouped <- gather(np.pah2.avg.spread2, variable, value, -Date)

np.pah2.sd.spread <- as.data.frame(spread(np.pah2.grouped.sd, variable, SD))
np.pah2.sd.spread2 <- np.pah2.sd.spread[, -1]
row.names(np.pah2.sd.spread2) <- np.pah2.sd.spread$Date
np.pah2.sd.spread2[is.na(np.pah2.sd.spread2)] <- 0
np.pah2.sd.spread2 <- mutate_all(np.pah2.sd.spread2, function(x) as.numeric(as.character(x)))
np.pah2.sd.spread2$Total_SD <- rowSums(np.pah2.sd.spread2, na.rm = T)
np.pah2.sd.spread2$Date <- row.names(np.pah2.sd.spread2)
np.pah2.sd.grouped <- gather(np.pah2.sd.spread2, variable, sd, -Date)

np.pah2.avg.grouped2 <- cbind(np.pah2.avg.grouped, np.pah2.sd.grouped[, -c(1,2)])
colnames(np.pah2.avg.grouped2) <- c("Date", "variable", "value", "sd")

SCOP.pah2.plot2 <- ggplot(np.pah2.avg.grouped2, aes(x = as.Date(Date), y = value)) +  
  geom_line(aes(color = variable, group = variable), linewidth = 0.5) + 
  geom_point(aes(color = variable), size = 3) + 
  geom_errorbar(aes(ymin = value-sd, ymax = value+sd, color = variable), width = 1.5) + 
  theme_classic() + 
  labs(x = "Date", y = "Avg. Concentration (ppb)") + 
  scale_x_date(breaks = as.Date(np.pah2.avg.grouped2$Date), date_labels =  "%Y-%m-%d") + 
  expand_limits(y = c(0, 125)) + 
  scale_y_continuous(breaks = seq(0, 125, 25)) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans", size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, colour = c("grey30", "grey30", "seagreen4", "firebrick3", rep("grey30", 5), "firebrick3"))) + 
  scale_color_manual(name = "Compound:", breaks = c("acenaphthene", "fluorene", "anthracene", "phenanthrene", "pyrene", "benz.a.anthracene", "Chrysene", "Total_PAH"), labels = c("Acenaphthene (LMW)", "Fluorene (LMW)", "Anthracene (LMW)", "Phenanthrene (LMW)", "Pyrene (HMW)", "Benz[a]anthracene (HMW)", "Chrysene (HMW)", "Total PAH Concentration"), values = c("darkgoldenrod3", "coral3", "darkslategray3", "darkslategray4", "grey80", "grey35", "black", "mediumorchid3"))
SCOP.pah2.plot2


####Export Figures####
setwd("G:/Figures/Environmental and Geochemical/")
getwd()

svg("fig1.svg", width = 25, height = 18)
png("fig1.png", width = 25, height = 18, units = "in", res = 600)
ggarrange(SCOP.temp.plot, SCOP.nitrate.plot, SCOP.POC.plot, SCOP.sal.plot, SCOP.SRP.plot, SCOP.chl.plot1, SCOP.CUTI33N.plot, SCOP.mwd.plot, ncol = 3, nrow = 3, align = "hv", common.legend = T, legend = "right", labels = c("a", "b", "c", "d", "e", "f", "g", "h"), font.label = list(size = 24, color = "black", face = "bold", family = "sans"))
dev.off()

png("fig2b-rev.png", width = 9, height = 6, units = "in", res = 600)
SCOP.pah2.plot2
dev.off()