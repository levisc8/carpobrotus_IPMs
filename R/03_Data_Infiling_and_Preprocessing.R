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

# Base graphics figures

# par(mfrow = c(2, 2))
# 
# hist(log(HavatseletRamets$Size),
#      main = '',
#      xlab = expression('Log ramet size distribution m' ^2),
#      breaks = 70)
# 
# hist(log(HavatseletGenets$Size),
#      main = '',
#      xlab = expression('Log genet size distribution m' ^2),
#      breaks = 10)
# 
# plot(Flower_N ~ Size, data = HavatseletRamets,
#      ylab = '# of Flowers per Ramet')
# 
# 
# plot(Flower_N ~ Size, data = HavatseletGenets,
#      ylab = '# of Flowers per Genet')

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
                     limit = c(-2, 7))

RametFecPlot <- ggplot(HavatseletRamets, 
                       aes(x = LogSize,
                           y = Flower_N)) + 
  geom_point(color = 'red') + 
  theme.bl + 
  scale_y_continuous('# of Flowers') + 
  scale_x_continuous(('Size of Ramets')) 

GenetFecPlot <- ggplot(HavatseletGenets,
                       aes(x = LogSize,
                           y = Flower_N)) + 
  geom_point(color = 'blue') + 
  theme.bl + 
  scale_y_continuous('# of Flowers') + 
  scale_x_continuous('Size of Genets',
                     limits = c(-2, 7))



pdf('Figures/Israel/Havatselet/Preliminary_Plots.pdf',
    width = 8,
    height = 8)
  grid.arrange(RametHist, GenetHist,
               RametFecPlot, GenetFecPlot,
               nrow = 2, ncol = 2)
dev.off()

png('Figures/Israel/Havatselet/Preliminary_Plots.png',
    width = 8,
    height = 8,
    units = 'in',
    res = 72)
  grid.arrange(RametHist, GenetHist,
               RametFecPlot, GenetFecPlot,
               nrow = 2, ncol = 2)
dev.off()

# If necessary, copy over png file so the collaboration outline can use it
# fs::file_copy('Figures/Israel/Havatselet/Preliminary_Plots.png',
#               '../Stellenbosch_Phys_Collab/R/Preliminary_Plots.png')
