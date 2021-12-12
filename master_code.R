#====================================
#SOUND ANALYSIS OF ORGAN PIPES
#====================================
#Lin Yangchen
#December 2021
#for personal non-commercial use only

rm(list = ls())



#====================================
#user settings
#====================================

#desired sample rate (downsampled if applicable)
downsamprate <- 96000

#where to extract stationary tone from recording
#in seconds
timerange <- c(4, 8)

#cut off power spectra at which harmonic
hcutoff <- 26

#frequency cutoff for stacked spectra
stackfreqcut <- 15000

#frequency band (Hz on each side) for flutter detection
flutterband <- c(50,50)

#no. of top notes to exclude from tuning and voicing analysis
keycut <- 47

#where to extract attack transient from recording
#in seconds
attacktime <- c(0.6, 2)

#cut off attack spectrograms at which harmonic
ahcutoff <- 11

#Google Sheet url of pipe data
dataurl <- 'docs.google.com/spreadsheets/d/1K_NoHORpwFQEOORFi61MO-eH3Msx34ED_nw8zD5yqhg#gid=0'

#====================================
#load libraries
#====================================


install.packages(setdiff(c(
'tuneR'
, 'stringr'
, 'lattice'
, 'latticeExtra'
, 'fBasics'
, 'extrafont'
, 'gsheet'
)
, rownames(installed.packages())))

require(tuneR)
require(stringr)
require(lattice)
require(latticeExtra)
require(fBasics)

require(extrafont)
#font_import()
#loadfonts(device = 'pdf')
#fonts()
font <- 'Courier New'

#plot settings
mgp <- c(2,0.5,0)
tck <- 0.01



#====================================
#import morphometric and metallurgical data
#====================================

require(gsheet)
pipedata <- read.csv(text = gsheet2text(dataurl, format = 'csv'))

#====================================
#run analysis on each replicate folder (set) of recordings
#====================================

folders <- list.dirs(recursive = F)

for(folder in folders)
    {
        file.copy('analysis_code.R', folder, overwrite = T)
        setwd(folder)
        source('analysis_code.R')
        setwd('..')
    }

