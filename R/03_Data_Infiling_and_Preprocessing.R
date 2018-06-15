# Basic shape file analysis

# Use sf for convenience in applying dplyr verbs to spatial objects.
# Notice that when aggregations are made, spatial information
# is still retained.

library(sf)
library(stringr)
library(ggplot2)
library(dplyr)
library(gridExtra)


infile_data('Polygons/Israel/Havatselet_Plants.shp',
            'HavatseletRamets',
            log.breaks = 60)

HavatseletGenets <- ramet_to_genet(HavatseletRamets,
                                   Clone_Of,
                                   size.breaks = 10)

# Generate info on reproductive status
HavatseletRamets$Repro <- ifelse(is.na(HavatseletRamets$Flower_N), 0, 1)

for(i in seq_len(dim(HavatseletGenets)[1])) {
  HavatseletGenets$Repro[i] <- ifelse(HavatseletGenets$Flower_N[i] == 0, 0, 1)
}

# First models. Fecundity!
RametFlowerGLM <- glm(Flower_N ~ LogSize, 
                      data = HavatseletRamets,
                      family = poisson())

GenetFlowerGLM <- glm(Flower_N ~ LogSize, 
                      data = HavatseletGenets,
                      family = poisson())

RametReproGLM <- glm(Repro ~ LogSize,
                     data = HavatseletRamets,
                     family = binomial())

GenetReproGLM <- glm(Repro ~ LogSize,
                     data = HavatseletGenets,
                     family = binomial())

message('Ramet flower #\n') # overdispersed
summary(RametFlowerGLM)
message('\nGenet flower #\n') # overdispersed
summary(GenetFlowerGLM)

message('\nRamet pr(Repro)\n') # not overdispersed!
summary(RametReproGLM)
message('\nGenet pr(Repro)\n') # overdispersed
summary(GenetFlowerGLM)

HavatseletGenets$PredPrRepro <- predict(GenetReproGLM,
                                        data.frame(LogSize = HavatseletGenets$LogSize),
                                        type = 'response')

HavatseletRamets$PredPrRepro <- predict(RametReproGLM,
                                        data.frame(LogSize = HavatseletRamets$LogSize),
                                        type = 'response')

HavatseletGenets$PredFlowers <- predict(GenetFlowerGLM,
                                        data.frame(LogSize = HavatseletGenets$LogSize),
                                        type = 'response')

HavatseletRamets$PredFlowers <- predict(RametFlowerGLM,
                                        data.frame(LogSize = HavatseletRamets$LogSize),
                                        type = 'response')

# ggplot2 + gridExtra

theme.bl <- theme(panel.background = element_rect(fill = NA,
                                                  color = 'black',
                                                  size = 1.25),
                  panel.grid = element_blank(),
                  axis.title.y = element_text(size = 14,
                                              margin = margin(t = 0,
                                                              l = 5,
                                                              r = 10, 
                                                              b = 0)),
                  axis.title.x = element_text(size = 14,
                                              margin = margin(t = 10, 
                                                              b = 5,
                                                              l = 0, 
                                                              r = 0)))
                  
RametHist <- ggplot(HavatseletRamets, aes(x = CleanBin)) + 
  geom_bar(width = 0.05,
           fill = 'red') + 
  theme.bl + 
  scale_y_continuous('# of Ramets') + 
  scale_x_continuous('Size of Ramets')

GenetHist <- ggplot(HavatseletGenets, aes(x = CleanBin)) + 
  geom_bar(width = 0.2,
           fill = 'blue') + 
  theme.bl + 
  scale_y_continuous('# of Genets') + 
  scale_x_continuous('Size of Genets',
                     limits = c(-2, 7))

RametFecPlot <- ggplot(HavatseletRamets, 
                       aes(x = LogSize,
                           y = Flower_N)) + 
  geom_point(color = 'red') + 
  geom_line(aes(y = PredFlowers),
            color = 'red') +
  theme.bl + 
  scale_y_continuous('# of Flowers') + 
  scale_x_continuous(('Size of Ramets')) 

GenetFecPlot <- ggplot(HavatseletGenets,
                       aes(x = LogSize,
                           y = Flower_N)) + 
  geom_point(color = 'blue') + 
  geom_line(aes(y = PredFlowers),
            color = 'blue') +
  theme.bl + 
  scale_y_continuous('# of Flowers') + 
  scale_x_continuous('Size of Genets',
                     limits = c(-2, 7))

RametReproPlot <- ggplot(HavatseletRamets,
                         aes(x = LogSize,
                             y = Repro)) + 
  geom_point(color = 'red',
             size = 1.25) + 
  geom_line(aes(y = PredPrRepro),
            color = 'red',
            linetype = 'dashed',
            size = 1.1,
            alpha = 0.5) + 
  theme.bl + 
  scale_y_continuous('Pr(Reproductive)') + 
  scale_x_continuous('Size of Ramets')



GenetReproPlot <- ggplot(HavatseletGenets,
                         aes(x = LogSize,
                             y = Repro)) + 
  geom_point(color = 'blue',
             size = 1.25) + 
  geom_line(aes(y = PredPrRepro),
            color = 'blue',
            linetype = 'dashed',
            size = 1.1,
            alpha = 0.5) + 
  theme.bl + 
  scale_y_continuous('Pr(Reproductive)') + 
  scale_x_continuous('Size of Genets',
                     limits = c(-2, 7))

pdf('Figures/Israel/Havatselet/Preliminary_Plots.pdf',
    width = 8,
    height = 8)
  grid.arrange(RametHist, GenetHist,
               RametFecPlot, GenetFecPlot,
               RametReproPlot, GenetReproPlot,
               nrow = 3, ncol = 2)
dev.off()

png('Figures/Israel/Havatselet/Preliminary_Plots.png',
    width = 8,
    height = 8,
    units = 'in',
    res = 72)
  grid.arrange(RametHist, GenetHist,
               RametFecPlot, GenetFecPlot,
               RametReproPlot, GenetReproPlot,
               nrow = 3, ncol = 2)
dev.off()

# If necessary, copy over png file so the collaboration outline can use it

if(fs::file_exists('../Stellenbosch_Phys_Collab/R/Preliminary_Plots.png')) {
  fs::file_delete('../Stellenbosch_Phys_Collab/R/Preliminary_Plots.png')
}
fs::file_copy('Figures/Israel/Havatselet/Preliminary_Plots.png',
              '../Stellenbosch_Phys_Collab/R/Preliminary_Plots.png')
