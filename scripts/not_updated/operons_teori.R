############### Plot af reference operons #########

#### Indl?s pakker ####

#library(ape)
#library(RFLPtools)
library(stringr)
library(tidyverse)
library(plyr)
library(dplyr)
#library(reshape2)
library(gplots)
#library(seqinr)
library(writexl)
#library(varhandle)
require('Gviz')
library(data.table)
require(stats)
#library(gridExtra)
library("readxl")
library(RColorBrewer)
library(xlsx)

##### WD #####
setwd("~/Desktop/AAU/8. semester/Projekt/Operon_plots")
#### Indl?s data ####

path_excel_file = "Query_figur.xlsx"
sheet_names = excel_sheets(path = path_excel_file)
list_all <- lapply(sheet_names, function(x) read_excel(path = path_excel_file, sheet = x, skip = 1))
list_all <- list_all[-20]
list_all_df = rbind.fill(list_all) 

list_all_df$Start = as.numeric(list_all_df$Start)
list_all_df$End = as.numeric(list_all_df$End)

#### Color after function ####

brewer.pal(n = 8, name = 'Oranges')


operons_combined = as_tibble(list_all_df)
operons_combined = operons_combined %>% 
  mutate(Color = case_when(
    ##### something with Transport / export  grå
    startsWith(Function, 'ABC') ~ '#dccede',
    startsWith(Function, 'PE') ~ '#C8DBC8',
    startsWith(Function, 'Export') ~ '#EDE9EB',
    ###### Modifications ## Greens
    startsWith(Function, 'Branch') ~ '#ffcb85',
    startsWith(Function, 'MOD') ~ '#deb887',
    ###### Break down #### 
    startsWith(Function, 'Hydrolase') ~ '#fcdbbb', 
    startsWith(Function, 'Lyase') ~ '#FBDADA',  
    ##### Synthases / GT  
    startsWith(Function, 'GT') ~ '#82ACB9',
    startsWith(Function, 'Polymerization') ~ '#a1c2a1',
    startsWith(Function, 'Prim_GT') ~ '#ECDED5',
    startsWith(Function, 'SY') ~ '#B3CFD0',
    
    ##### Other   
    startsWith(Function, 'NA') ~ '#FFFFFF', #Grey
    startsWith(Function, 'PS') ~ '#D68A8A',  #Purple
    startsWith(Function, 'REG') ~ '#f5f5c6'  # yellow
  ))


operons_combined = operons_combined %>% 
  mutate(Correct_names = case_when(
    ##### something with Transport / export  BLUES
    startsWith(Function, 'ABC') ~ 'ABC Transporter',
    startsWith(Function, 'PE') ~ 'Polymerisation and Export',
    startsWith(Function, 'Export') ~ 'Export',
    ###### Modifications ## Greens
    startsWith(Function, 'Branch') ~ 'Branching',
    startsWith(Function, 'Hydrolase') ~ 'Hydrolase',
    startsWith(Function, 'Lyase') ~ 'Lyase',
    startsWith(Function, 'MOD') ~ 'Modifications',
    ##### Synthases / GT  Oranges
    startsWith(Function, 'GT') ~ 'GT',
    startsWith(Function, 'Polymerization') ~ 'Polymerisation',
    startsWith(Function, 'SY') ~ 'Synthase',
    startsWith(Function, 'Prim_GT') ~ 'Priming GT',
    startsWith(Function, 'PCP') ~ 'PCP',
    ##### Other   
    startsWith(Function, 'NA') ~ 'Unknown', #Grey
    startsWith(Function, 'PS') ~ 'Precurser Synthase',  #Purple
    startsWith(Function, 'REG') ~ 'Regulation'  # yellow
  ))




#Structure data from operons combined to GRanges format ####

#plot_parameters
p_width = 10
p_showFeatureId = T  
p_min.height = 10 
p_cex=0.4
p_lwd = 1
p_arrowHeadWidth=30
p_fontcolor.item="grey1"
p_shape = 'arrow'
p_just.group = 'above'
p_rotation.item=0
p_fontface.group = 12
p_fontfamily.group ='mono'
p_fontcolor.group ="grey1"


#### data_track --------------------------------------------------------------



#Alginate
data_Alginate = filter(operons_combined, Polysaccheride == 'Alginate')
GR_data_Alginate <- with(data_Alginate, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Alginate <- AnnotationTrack(GR_data_Alginate, name = "Alginate",  
                                       width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                       cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                       fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                       stacking = 'dense', id = data_Alginate$Genename, fill = data_Alginate$Color)

#Cellulose1
data_Cellulose1 = filter(operons_combined, Polysaccheride == 'Cellulose1')
GR_data_Cellulose1 <- with(data_Cellulose1, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Cellulose1 <- AnnotationTrack(GR_data_Cellulose1, name = "Cellu- lose 1",  
                                         width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                         cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35', 
                                         fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                         stacking = 'dense', id = data_Cellulose1$Genename, fill = data_Cellulose1$Color)

#Cellulose2
data_Cellulose2 = filter(operons_combined, Polysaccheride == 'Cellulose2')
GR_data_Cellulose2 <- with(data_Cellulose2, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Cellulose2 <- AnnotationTrack(GR_data_Cellulose2, name = "Cellu- lose 2", 
                                         width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                         cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                         fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                         stacking = 'dense', id = data_Cellulose2$Genename, fill = data_Cellulose2$Color)

#Curdlan
data_Curdlan = filter(operons_combined, Polysaccheride == 'Curdlan')
GR_data_Curdlan <- with(data_Curdlan, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Curdlan <- AnnotationTrack(GR_data_Curdlan, name = "Curdlan", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_Curdlan$Genename, fill = data_Curdlan$Color)

#Diutan
data_Diutan = filter(operons_combined, Polysaccheride == 'Diutan')
GR_data_Diutan <- with(data_Diutan, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Diutan <- AnnotationTrack(GR_data_Diutan, name = "Diutan",  
                                     width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                     cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                     fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                     stacking = 'dense', id = data_Diutan$Genename, fill = data_Diutan$Color)

#Gellan
data_Gellan = filter(operons_combined, Polysaccheride == 'Gellan')
GR_data_Gellan <- with(data_Gellan, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Gellan <- AnnotationTrack(GR_data_Gellan, name = "Gellan", 
                                     width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                     cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                     fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                     stacking = 'dense', id = data_Gellan$Genename, fill = data_Gellan$Color)

#HAmulto
data_HAmulto = filter(operons_combined, Polysaccheride == 'HAmulto')
GR_data_HAmulto <- with(data_HAmulto, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_HAmulto <- AnnotationTrack(GR_data_HAmulto, name = "HA 1", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_HAmulto$Genename, fill = data_HAmulto$Color)

#HApyogenes
data_HApyogenes = filter(operons_combined, Polysaccheride == 'HApyogenes')
GR_data_HApyogenes <- with(data_HApyogenes, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_HApyogenes <- AnnotationTrack(GR_data_HApyogenes, name = "HA 3", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_HApyogenes$Genename, fill = data_HApyogenes$Color)

#HAequi
data_HAequi = filter(operons_combined, Polysaccheride == 'HAequi')
GR_data_HAequi <- with(data_HAequi, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_HAequi <- AnnotationTrack(GR_data_HAequi, name = "HA 2", 
                                     width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                     cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                     fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                     stacking = 'dense', id = data_HAequi$Genename, fill = data_HAequi$Color)

#Pel
data_Pel = filter(operons_combined, Polysaccheride == 'Pel')
GR_data_Pel <- with(data_Pel, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Pel <- AnnotationTrack(GR_data_Pel, name = "Pel", 
                                  width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                  cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                  fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                  stacking = 'dense', id = data_Pel$Genename, fill = data_Pel$Color)

#psl
data_psl = filter(operons_combined, Polysaccheride == 'psl')
GR_data_psl <- with(data_psl, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_psl <- AnnotationTrack(GR_data_psl, name = "psl", 
                                  width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                  cex = p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                  fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                  stacking = 'dense', id = data_psl$Genename, fill = data_psl$Color)

#S88
data_S88 = filter(operons_combined, Polysaccheride == 'S88')
GR_data_S88 <- with(data_S88, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_S88 <- AnnotationTrack(GR_data_S88, name = "S88", 
                                  width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                  cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                  fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                  stacking = 'dense', id = data_S88$Genename, fill = data_S88$Color)

#Salecan
data_Salecan = filter(operons_combined, Polysaccheride == 'Salecan')
GR_data_Salecan <- with(data_Salecan, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Salecan <- AnnotationTrack(GR_data_Salecan, name = "Salecan", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_Salecan$Genename, fill = data_Salecan$Color)


#Succinoglycan
data_Succinoglycan = filter(operons_combined, Polysaccheride == 'Succinoglycan')
GR_data_Succinoglycan <- with(data_Succinoglycan, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Succinoglycan <- AnnotationTrack(GR_data_Succinoglycan, name = "Succino- glycan", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_Succinoglycan$Genename, fill = data_Succinoglycan$Color)


#Xanthan
data_Xanthan = filter(operons_combined, Polysaccheride == 'Xanthan')
GR_data_Xanthan <- with(data_Xanthan, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_Xanthan <- AnnotationTrack(GR_data_Xanthan, name = "Xanthan", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_Xanthan$Genename, fill = data_Xanthan$Color)


#PNAG_pga
data_PNAG_pga = filter(operons_combined, Polysaccheride == 'PNAG_pga')
GR_data_PNAG_pga <- with(data_PNAG_pga, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_PNAG_pga <- AnnotationTrack(GR_data_PNAG_pga, name = "PNAG (pga)", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_PNAG_pga$Genename, fill = data_PNAG_pga$Color)



#PNAG_ica
data_PNAG_ica = filter(operons_combined, Polysaccheride == 'PNAG_ica')
GR_data_PNAG_ica <- with(data_PNAG_ica, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_PNAG_ica <- AnnotationTrack(GR_data_PNAG_ica, name = "PNAG (ica)", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_PNAG_ica$Genename, fill = data_PNAG_ica$Color)


#PNAG_eps
data_PNAG_eps = filter(operons_combined, Polysaccheride == 'PNAG_eps')
GR_data_PNAG_eps <- with(data_PNAG_eps, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_PNAG_eps <- AnnotationTrack(GR_data_PNAG_eps, name = "PNAG (eps)", 
                                      width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                      cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35',
                                      fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                      stacking = 'dense', id = data_PNAG_eps$Genename, fill = data_PNAG_eps$Color)



#NulO
data_NulO = filter(operons_combined, Polysaccheride == 'NulO')
GR_data_NulO <- with(data_NulO, GRanges(Operon, IRanges(Start, End), Strand, id=Genename, Color = Color))
data_track_NulO <- AnnotationTrack(GR_data_NulO, name = "NulO", 
                                       width = p_width, showFeatureId = p_showFeatureId, min.height= p_min.height, p_fontcolor.group ="grey1",
                                       cex= p_cex, arrowHeadMaxWidth= p_arrowHeadWidth, fontface.group = p_fontface.group, fontfamily = p_fontfamily.group, background.title = 'grey35', 
                                       fontcolor.item=p_fontcolor.item, shape = p_shape, just.group = p_just.group, lwd = p_lwd, rotation.item = p_rotation.item,
                                       stacking = 'dense', id = data_NulO$Genename, fill = data_NulO$Color)










#### plotting ----------------------------------------------------------------

## Create legend


operons_combined_1 = as_tibble(operons_combined)
operons_combined_1 = filter(operons_combined_1, Polysaccheride == 'Alginate' | Polysaccheride == 'Cellulose1' | Polysaccheride == 'Cellulose2' |
                                Polysaccheride == 'Curdlan' | Polysaccheride == 'HAequi' | Polysaccheride == 'HApyogenes' | 
                                Polysaccheride == 'HAmulto' | Polysaccheride == 'Pel' | Polysaccheride == 'NulO' | Polysaccheride == 'HAequi' | 
                                Polysaccheride == 'PNAG_ica' | Polysaccheride == 'PNAG_pga' | Polysaccheride == 'PNAG_eps')


legend_1 <- operons_combined_1 %>% 
  distinct(Function, .keep_all = TRUE) %>%
  select(Function, Color, Correct_names) %>%
  arrange(Correct_names)


operons_combined_23 = as_tibble(operons_combined)
operons_combined_23 = filter(operons_combined_23, Polysaccheride == 'S88' | Polysaccheride == 'Gellan' | Polysaccheride == 'Diutan' |
                               Polysaccheride == 'psl' | Polysaccheride == 'Salecan' | Polysaccheride == 'Succinoglycan' | 
                               Polysaccheride == 'Xanthan')


legend_23 <- operons_combined_23 %>% 
  distinct(Function, .keep_all = TRUE) %>%
  select(Function, Color, Correct_names) %>%
  arrange(Correct_names)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", legend = pull(legend_1, Correct_names) , col = pull(legend_1, Color),
       pch=15, pt.cex=2, cex=0.5, bty='n', text.col= 'grey1', ncol = 5, 
       x.intersp = 1.2, y.intersp = 2)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", legend = pull(legend_23, Correct_names), col = pull(legend_23, Color),
       pch=15, pt.cex=2, cex=0.5, bty='n', text.col= 'grey1', ncol = 7, x.intersp = 1.2, y.intersp = 2)


######### Short operons ###############

####Make bounderies and bp line in plot ####

ref_2 <- GRanges('Operon', IRanges(0, 17000))
ref_track_2 <- GenomeAxisTrack(ref_2, lwd=1, fontsize=10, labelPos="above", 
                               distFromAxis=2, ticksAt= c(seq.int(2000, 18000, 2000)))

plot_operon_ref_2 = plotTracks(c(ref_track_2, 
             data_track_Alginate, data_track_Cellulose1, data_track_Cellulose2, 
             data_track_Curdlan, data_track_HAmulto, data_track_HApyogenes, data_track_HAequi, data_track_Pel, data_track_NulO,
             data_track_PNAG_pga, data_track_PNAG_ica, data_track_PNAG_eps),
              cex.main = 0.5)

##### Large operons ##################

ref_1 <- GRanges('Operon', IRanges(0, 30000))
ref_track_1 <- GenomeAxisTrack(ref_1, lwd=1, fontsize=10, 
                               labelPos="above", ticksAt= c(seq.int(2000, 30000, 4000), cex.title=2))

plot_operon_ref_1 = plotTracks(c(ref_track_1, 
                                 data_track_Diutan, data_track_Gellan, 
                                 data_track_S88, data_track_Succinoglycan), 
                                cex.main = 0.5)

########## Small NON synthase

ref_3 <- GRanges('Operon', IRanges(0, 16000))
ref_track_3 <- GenomeAxisTrack(ref_3, lwd=1, fontsize=10, 
                               labelPos="above", ticksAt= c(seq.int(2000, 18000, 2000), cex.title=2))

plot_operon_ref_3 = plotTracks(c(ref_track_3, 
                                 data_track_psl, data_track_Salecan, data_track_Xanthan), 
                               cex.main = 0.5)


###


