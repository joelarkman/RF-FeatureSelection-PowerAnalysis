#################################################################
######################## Data Cleaning ##########################
#################################################################

# Load packages to access GEO database
library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE62932", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "33224434223431434111241441323323112243243233224421234433131331220000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- factor(sml,labels = c("Control","Stage_1","Stage_2","Stage_3","Stage_4"))


# Take known metadata from GEO database
raw_metadata <- data.frame('GEO_ID'=colnames(ex),'Metadata'=fl)
# Import metadata file.
match_metadata <- read.csv('key.csv')
# Match database and uploaded metadata.
index <-match(raw_metadata$GEO_ID, match_metadata$GEO_ID)
# Create final metadata dataframe. 
metadata <- data.frame(match_metadata[index,], 'Metadata'=raw_metadata$Metadata)

# Load actual genomic data.
matched_data <- read.csv('GSE62932_microarray_matched.csv', row.names = 1)

# Match samples to metadata.
index <- match(paste('X',metadata$Matched_ID, sep = ''),colnames(matched_data))

# Store final matched dataset.
matched_data <- t(matched_data[,index])


#################################################################
##################### Classification plot #######################
#################################################################

# Extract stage 1 CRC data.
xdata <- matched_data[metadata$Metadata %in% c('Control','Stage_1'),]
ydata <- droplevels(as.factor(metadata$Metadata[metadata$Metadata %in% c('Control','Stage_1')]))
stage_1 <- data.frame(y=ydata,xdata)

# Extract stage 2 CRC data.
xdata <- matched_data[metadata$Metadata %in% c('Control','Stage_2'),]
ydata <- droplevels(as.factor(metadata$Metadata[metadata$Metadata %in% c('Control','Stage_2')]))
stage_2 <- data.frame(y=ydata,xdata)

# Extract stage 3 CRC data.
xdata <- matched_data[metadata$Metadata %in% c('Control','Stage_3'),]
ydata <- droplevels(as.factor(metadata$Metadata[metadata$Metadata %in% c('Control','Stage_3')]))
stage_3 <- data.frame(y=ydata,xdata)

# Extract stage 4 CRC data.
xdata <- matched_data[metadata$Metadata %in% c('Control','Stage_4'),]
ydata <- droplevels(as.factor(metadata$Metadata[metadata$Metadata %in% c('Control','Stage_4')]))
stage_4 <- data.frame(y=ydata,xdata)


# Store names of features chosen for stage 1.
variables_S1 <- c('CXCR7\n(Effect: 4.32, Group size: 4)',
                  'PCNA\n(Effect: 3.24, Group size: 3)',
                  'ABP1\n(Effect: 3.06, Group size: 12)')

# Extract data for feartures chosen for stage 1.
s1_v1 <- data.frame(y=stage_1$y,x=stage_1$CXCR7, v=variables_S1[1], stage=1)
s1_v2 <- data.frame(y=stage_1$y,x=stage_1$PCNA, v=variables_S1[2], stage=1)
s1_v3 <- data.frame(y=stage_1$y,x=stage_1$ABP1, v=variables_S1[3], stage=1)

# Store names of features chosen for stage 2.
variables_S2 <- c('PRV1\n(Effect: 4.19, Group size: 7)',
                  'MYOT\n(Effect: 3.08, Group size: 8)',
                  'CRABP1\n(Effect: 2.39, Group size: 1)',
                  'CLDN1\n(Effect: 2.21, Group size: 5)',
                  'BMP2\n(Effect: 1.98, Group size: 3)',
                  'SUMO2\n(Effect: 1.91, Group size: 1)')

# Extract data for feartures chosen for stage 2.
s2_v1 <- data.frame(y=stage_2$y,x=stage_2$PRV1, v=variables_S2[1], stage=2)
s2_v2 <- data.frame(y=stage_2$y,x=stage_2$MYOT, v=variables_S2[2], stage=2)
s2_v3 <- data.frame(y=stage_2$y,x=stage_2$CRABP1, v=variables_S2[3], stage=2)
s2_v4 <- data.frame(y=stage_2$y,x=stage_2$CLDN1, v=variables_S2[4], stage=2)
s2_v5 <- data.frame(y=stage_2$y,x=stage_2$BMP2, v=variables_S2[5], stage=2)
s2_v6 <- data.frame(y=stage_2$y,x=stage_2$SUMO2, v=variables_S2[6], stage=2)

# Store names of features chosen for stage 3.
variables_S3 <- c('PRV1\n(Effect: 4.78, Group size: 4)',
                  'HMP19\n(Effect: 2.99, Group size: 5)',
                  'HMMR\n(Effect: 2.94, Group size: 4)',
                  'HNRPA3\n(Effect: 2.9, Group size: 3)',
                  'HNRNPM\n(Effect: 2.73, Group size: 3)',
                  'TIMP1\n(Effect: 2.48, Group size: 7)',
                  'CRABP1\n(Effect: 2.47, Group size: 2)',
                  'CXCL16\n(Effect: 2.22, Group size: 2)',
                  'RFWD3\n(Effect: 2.08, Group size: 3)',
                  'IFITM1\n(Effect: 1.99, Group size: 1)',
                  'PPIB\n(Effect: 1.71, Group size: 2)',
                  'HRASLS3\n(Effect: 1.65, Group size: 1)')

# Extract data for feartures chosen for stage 3.
s3_v1 <- data.frame(y=stage_3$y,x=stage_3$PRV1, v=variables_S3[1], stage=3)
s3_v2 <- data.frame(y=stage_3$y,x=stage_3$HMP19, v=variables_S3[2], stage=3)
s3_v3 <- data.frame(y=stage_3$y,x=stage_3$HMMR, v=variables_S3[3], stage=3)
s3_v4 <- data.frame(y=stage_3$y,x=stage_3$HNRPA3, v=variables_S3[4], stage=3)
s3_v5 <- data.frame(y=stage_3$y,x=stage_3$HNRNPM, v=variables_S3[5], stage=3)
s3_v6 <- data.frame(y=stage_3$y,x=stage_3$TIMP1, v=variables_S3[6], stage=3)
s3_v7 <- data.frame(y=stage_3$y,x=stage_3$CRABP1, v=variables_S3[7], stage=3)
s3_v8 <- data.frame(y=stage_3$y,x=stage_3$CXCL16, v=variables_S3[8], stage=3)
s3_v9 <- data.frame(y=stage_3$y,x=stage_3$RFWD3, v=variables_S3[9], stage=3)
s3_v10 <- data.frame(y=stage_3$y,x=stage_3$IFITM1, v=variables_S3[10], stage=3)
s3_v11 <- data.frame(y=stage_3$y,x=stage_3$PPIB, v=variables_S3[11], stage=3)
s3_v12 <- data.frame(y=stage_3$y,x=stage_3$HRASLS3, v=variables_S3[12], stage=3)

# Store names of features chosen for stage 4.
variables_S4 <- c('PRV1\n(Effect: 6.1, Group size: 4)',
                  'MYOT\n(Effect: 4.56, Group size: 4)',
                  'S100A11\n(Effect: 2.99, Group size: 5)',
                  'C10orf22\n(Effect: 2.65, Group size: 8)',
                  'FTH1\n(Effect: 2.56, Group size: 3)',
                  'BSG\n(Effect: 1.91, Group size: 2)',
                  'CLDN1\n(Effect: 1.9, Group size: 2)',
                  'SLC4A4\n(Effect: 1.74, Group size: 2)',
                  'MMP13\n(Effect: 1.19, Group size: 1)')

# Extract data for feartures chosen for stage 4.
s4_v1 <- data.frame(y=stage_4$y,x=stage_4$PRV1, v=variables_S4[1], stage=4)
s4_v2 <- data.frame(y=stage_4$y,x=stage_4$MYOT, v=variables_S4[2], stage=4)
s4_v3 <- data.frame(y=stage_4$y,x=stage_4$S100A11, v=variables_S4[3], stage=4)
s4_v4 <- data.frame(y=stage_4$y,x=stage_4$C10orf22, v=variables_S4[4], stage=4)
s4_v5 <- data.frame(y=stage_4$y,x=stage_4$FTH1, v=variables_S4[5], stage=4)
s4_v6 <- data.frame(y=stage_4$y,x=stage_4$BSG, v=variables_S4[6], stage=4)
s4_v7 <- data.frame(y=stage_4$y,x=stage_4$CLDN1, v=variables_S4[7], stage=4)
s4_v8 <- data.frame(y=stage_4$y,x=stage_4$SLC4A4, v=variables_S4[8], stage=4)
s4_v9 <- data.frame(y=stage_4$y,x=stage_4$MMP13, v=variables_S4[9], stage=4)

# Combine all data.
data <- rbind(s1_v1,s1_v2,s1_v3,
              s2_v1,s2_v2,s2_v3,s2_v4,s2_v4,s2_v5,s2_v6,
              s3_v1,s3_v2,s3_v3,s3_v4,s3_v5,s3_v6,s3_v7,s3_v8,s3_v9,s3_v10,s3_v11,s3_v12,
              s4_v1,s4_v2,s4_v3,s4_v4,s4_v5,s4_v6,s4_v7,s4_v8,s4_v9)

# Load ggpubr package.
library(ggpubr)

# Create plot.
ggplot(data, aes(x=y, y=x, fill=y)) +
  geom_boxplot() +
  facet_wrap(~v, scales = 'free', nrow = 4) +
  theme_bw() +
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  ylab('Variable Values') +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = 1, label.x.npc = 1, size =3)
  scale_fill_manual(values = c('#DF6F6F','#9CD1DE'))

# Store plot.
pdf('stage4.pdf', width=2*9.2, height=2)
ggplot(data[which(data$stage==4),], aes(x=y, y=x, fill=y)) +
  geom_boxplot() +
  facet_wrap(~v, scales = 'free', nrow = 1) +
  theme_bw() +
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  ylab('Variable Values') +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = 1, label.x.npc = 1, size =3) +
  scale_fill_manual(values = c('#DF6F6F','#F0E882'))
dev.off()

#BE92E1 - purple 3
#F0E882 - yellow 4
#AFE192 - green 2
#9CD1DE - blue 1


