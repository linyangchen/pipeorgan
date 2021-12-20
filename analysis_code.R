#====================================
#SOUND ANALYSIS OF ORGAN PIPES
#====================================
#Lin Yangchen
#December 2021
#for personal non-commercial use only



#====================================
#import wav files
#====================================

#if stereo, takes left channel
#if mono right channel, copies data to left channel and uses left channel

print(noquote(paste('importing wav files from', folder, 'folder')))

#numerical file names correspond to unique pipe IDs in pipedata
soundfiles <- list.files(pattern = '\\.wav$', ignore.case = T)
nclips <- length(soundfiles)

clipid <- str_remove(soundfiles, ".wav")
clipid <- str_remove(clipid, ".WAV")
clipid <- as.numeric(clipid)


sounddata    <- list(NULL)
maxtimestep  <- 0
minamplitude <- 0
maxamplitude <- 0
samprate     <- NULL
bitdepth     <- NULL
for(i in 1:nclips)
    {
        rawdata <- readWave(soundfiles[i])
        
        if(i == 1 &
        rawdata@samp.rate > downsamprate)
            {print(noquote(paste(
            'downsampling from', rawdata@samp.rate/1000,
            'kHz to', downsamprate/1000, 'kHz'
            )))}
        
        print(noquote(paste('pipe', clipid[i])))
        
        sampledata <- downsample(rawdata, downsamprate)
        
        if(length(sampledata@left) == 0)
            {
                sampledata@left <- sampledata@right
                sampledata@right <- numeric(0)
            }

        sounddata[[clipid[i]]] <- sampledata

        
        #no. of samples in longest recording
        maxtimestep <- max(maxtimestep, length(sounddata[[clipid[i]]]@left))

        #range of amplitude recorded across clips
        minamplitude <- min(minamplitude, min(sounddata[[clipid[i]]]@left))
        maxamplitude <- max(maxamplitude, max(sounddata[[clipid[i]]]@left))
        
        samprate <- c(samprate, sounddata[[clipid[i]]]@samp.rate)
        if(length(samprate) > 1 & var(samprate)!=0)
            {stop('clips have different sample rates')}

        bitdepth <- c(bitdepth, sounddata[[clipid[i]]]@bit)
        if(length(bitdepth) > 1 & var(bitdepth)!=0)
            {stop('clips have different bit depths')}

    }


maxtimestep <- maxtimestep - 1
samprate <- samprate[1]
bitdepth <- bitdepth[1]

save(sounddata, file = 'sounddata.RData')




#====================================
#waveforms
#====================================

print(noquote(paste('plotting pipe waveforms in', folder, 'folder')))

xmax <- maxtimestep/samprate
pdf(height = 5, width = 5, 'waveforms.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

    for(i in clipid)
        {
            print(noquote(paste('pipe', i)))
            
            x <- 0:(length(sounddata[[i]]@left) - 1)
            x <- x/samprate
            
            plot(x, sounddata[[i]]@left/maxamplitude
            , type = 'l', lwd = 0.25
            , xlim = c(0, xmax), xaxs = 'i'
            , ylim = c(minamplitude/maxamplitude, 1)
            , mgp = mgp, tck = tck, axes = F
            , xlab = 'time (s)', ylab = 'normalized sound pressure'
            )
            
            axis(1, mgp = mgp, tck = tck)
            axis(2, mgp = mgp, tck = tck)
            box(bty = 'L')
            
            text(xmax, 1,
            paste('pipe', i,
            samprate/1000, 'kHz', bitdepth, 'bit'),
            pos = 2)
        }



dev.off()
embedFonts('waveforms.pdf')





#====================================
#stationary spectra
#====================================

print(noquote(paste(
'volume, fundamental frequencies and spectra/melfcc in',
folder, 'folder')))

samplerange <- timerange/xmax*maxtimestep

statspectra  <- list(NULL)
statcepstra  <- list(NULL)
avmelfcc     <- matrix(NA, nrow = nrow(pipedata), ncol = numcep - 1)
minspecamp   <- NULL
maxspecamp   <- NULL
fundamentals <- rep(NA, nrow(pipedata))
fundamps     <- rep(NA, nrow(pipedata))
rms          <- rep(NA, nrow(pipedata))
ampvar       <- rep(NA, nrow(pipedata))
ampdistrib   <- list(NULL)
scentroid    <- rep(NA, nrow(pipedata))

for(i in clipid)
    {
        print(noquote(paste('calculating pipe', i)))
        
        subdata <- sounddata[[i]]
        subdata@left <- subdata@left[samplerange[1]:samplerange[2]]
        
        rms[i] <- sqrt(mean(subdata@left^2))
        ampvar[i] <- var(subdata@left)
        ampdistrib[[i]] <- density(subdata@left)
        
        hs <- spectrum(subdata@left
        #, span = #smoothing
        , plot = F)
        
        statspectra[[i]] <- hs
        
        minspecamp <- min(minspecamp, min(hs$spec))
        maxspecamp <- max(maxspecamp, max(hs$spec))
        
        #need to calculate only once
        if(i == clipid[1])
        {
            freqseq <- seq(0, samprate/2, length.out = length(hs$freq))
            
            centind <- which(freqseq <= centfreqcut)
        }
        
        #determine fundamental frequency
        
        ref.freq <- pipedata$freq.eq.Hz[i]
        
        #look for peak
        freqinds <- which(freqseq > ref.freq*0.7 & freqseq < ref.freq*1.8)
        peak <- max(hs$spec[freqinds])
        fundamps[i] <- peak
        peakind <- which(hs$spec[freqinds] == peak)
        fundfreq <- freqseq[freqinds[peakind]]
        fundamentals[i] <- fundfreq
        
        #spectral centroid
        scentroid[i] <- weighted.mean(freqseq[centind],
        w = hs$spec[centind])
        
        #melfcc
        
        cepstrum <- melfcc(subdata
        , sr = samprate
        , maxfreq = 20000
        , numcep = numcep
        )
        
        statcepstra[[i]] <- cepstrum
        avmelfcc[i,] <- apply(cepstrum, 2, mean)[-1]
    }

ampvar <- ampvar/max(ampvar, na.rm = T)

save(statspectra, file = 'statspectra.RData')
save(statcepstra, file = 'statcepstra.RData')
save(avmelfcc, file = 'avmelfcc.RData')
save(fundamentals, file = 'fundamentals.RData')
save(fundamps, file = 'fundamps.RData')
save(rms, file = 'rms.RData')
save(ampvar, file = 'ampvar.RData')
save(ampdistrib, file = 'ampdistrib.RData')
save(scentroid, file = 'scentroid.RData')



#probability density distributions of sound pressure time series

pdf(height = 5, width = 5, 'time_pressure_distributions.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

for(i in clipid)
    {
        
        plot(ampdistrib[[i]]
        , zero.line = F
        , xaxs = 'i'
        , mgp = mgp, tck = tck, axes = F
        #, xlab = 'sound pressure', ylab = 'probability density'
        , xlab = '', ylab = ''
        , main = ''
        )
        
        mtext(paste('pipe', i))
        
        #axis(1, mgp = mgp, tck = -0.01)
        #axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')

    }


dev.off()
embedFonts('time_pressure_distributions.pdf')





#plot individual stationary spectra
pdf(height = 5, width = 5, 'stationary_spectra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1)
)


for(i in clipid)
    {
        print(noquote(paste('plotting pipe', i, 'spectrum')))
        
        #plot frequency cut-off at specified harmonic
        if(!is.na(hcutoff))
        {
            freqcutoff <- fundamentals[i]*hcutoff
        } else
        {
            freqcutoff <- samprate/2
        }
        freqcutind <- which(freqseq <= freqcutoff)
        
        ylim = c(log10(minspecamp), log10(maxspecamp))*10
        
        plot(freqseq[freqcutind], log10(statspectra[[i]]$spec)[freqcutind]*10
        , type = 'l', lwd = 0.5
        , xlim = c(0, freqcutoff)
        , ylim = ylim
        , xaxs = 'i'
        , mgp = mgp, tck = tck, axes = F
        , xlab = 'frequency (Hz)', ylab = 'dB'
        )
        
        axis(1, mgp = mgp, tck = tck)
        axis(2, mgp = mgp, tck = tck)
        box(bty = 'L')
        
        text(fundamentals[i],
        log10(statspectra[[i]]$spec)[which(freqseq == fundamentals[i])]*10,
        sprintf('%0.2f Hz', round(fundamentals[i], digits = 2)),
        #labels = paste(round(fundamentals[i], digits = 2), 'Hz'),
        pos = 4
        )
        
        text(freqcutoff, max(ylim),
        paste('pipe', i, 'to harmonic', hcutoff - 1)
        , pos = 2)

    }

dev.off()
embedFonts('stationary_spectra.pdf')




#plot individual stationary cepstra
pdf(height = 5, width = 5, 'stationary_cepstra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1)
)

for(i in clipid)
    {
        print(noquote(paste('plotting pipe', i, 'cepstrum')))

        image(statcepstra[[i]][,-1]
        , col = timPalette(10000)
        , xlab = 'time', ylab = 'melfcc'
        , mgp = mgp, tck = -0.01, axes = F
        )

        axis(1, mgp = mgp, tck = -0.01)
        box(bty = 'L')
        
        mtext(paste('pipe', i))

    }

dev.off()
embedFonts('stationary_cepstra.pdf')

#====================================
#flutter detection
#====================================

print(noquote(paste('detecting pipe flutter in',
folder, 'folder')))

pdf(height = 5, width = 5, 'flutter.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

for(i in clipid)
    {
        print(noquote(paste('pipe', i)))
        
        #freq cut-off below/above fundamental
        bandleft  <- fundamentals[i] - flutterband
        bandright <- fundamentals[i] + flutterband
        freqcutind <- which(freqseq > bandleft &
        freqseq < bandright)
        
        yvals <- log10(statspectra[[i]]$spec)[freqcutind]*10
        plot(freqseq[freqcutind], yvals
        , type = 'l', lwd = 0.5
        , xaxs = 'i'
        , mgp = mgp, tck = tck, axes = F
        , xlab = 'frequency (Hz)', ylab = 'dB'
        )
        
        axis(1, mgp = mgp, tck = tck)
        axis(2, mgp = mgp, tck = tck)
        box(bty = 'L')
        
        text(max(freqseq[freqcutind]), max(yvals),
        paste('pipe', i)
        , pos = 2)
    
    }
    
dev.off()
embedFonts('flutter.pdf')








#====================================
#voicing
#====================================


#3D stack of stationary power spectra of pipes in each stop
print(noquote(paste('plotting stacked voicing spectra in',
folder, 'folder')))

pdf(height = 5, width = 5, 'voicing_spectra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

stops <- unique(pipedata$stop[clipid])
spec3d <- list(NULL)
for(i in stops)
    {
        print(noquote(paste(i, 'spectra')))
        
        pipeid <- clipid[pipedata$stop[clipid] %in% i]
        
        #use only the specified no. of notes from the bottom of the keyboard
        keys <- pipedata$key[pipeid]
        keys <- keys[which(keys <= keycut)]

        #frequency cut-off at specified harmonic
        if(is.na(stackfreqcut))
            {
                stackfreqcut <- samprate/2
            }
        stackfreqcutind <- which(freqseq <= stackfreqcut)

        numfreqbin <- length(freqseq[stackfreqcutind])
        
        spec3d[[i]] <- matrix(NA, ncol = max(keys), nrow = numfreqbin)
        
        for(j in 1:length(keys))
            {
                #retrieve spectrum
                spec3d[[i]][,keys[j]] <- statspectra[[pipeid[j]]]$spec[stackfreqcutind]
            }
        
        #plot spectra
        
        image(freqseq[stackfreqcutind], 1:max(keys)
        , log10(spec3d[[i]])*10
        , col = timPalette(10000)
        , mgp = mgp, tck = -0.01, axes = F
        , xlab = 'frequency (Hz)', ylab = 'note'
        , main = i
        )
        
        axis(1, mgp = mgp, tck = -0.01)
        axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')

    }


dev.off()
embedFonts('voicing_spectra.pdf')



#3D stack of stationary melfcc averaged over time

print(noquote(paste('plotting stacked voicing cepstra in',
folder, 'folder')))

pdf(height = 5, width = 9, 'voicing_cepstra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1)
, mfrow = c(1, 2))
        
for(i in stops)
    {
        
        image(
        pipedata$key[which(pipedata$stop == i)]
        , 1:(numcep - 1)
        , avmelfcc[which(pipedata$stop == i),]
        , col = timPalette(10000)
        , mgp = mgp, tck = -0.01, axes = F
        , xlab = 'note', ylab = 'stationary melfcc'
        )
    
        axis(1, mgp = mgp, tck = -0.01)
        box(bty = 'L')
    
        mtext(i)    
        
    }


dev.off()
embedFonts('voicing_cepstra.pdf')



#====================================
#tuning and volume
#====================================

print(noquote(paste('analyzing temperament of each stop in',
folder, 'folder')))

#use only the specified no. of notes from the bottom of the keyboard
exclind <- NULL
for(i in stops)
    {
        stopind <- which(pipedata$stop == i)
        exclind <- c(exclind, stopind[1:keycut])
    }
truncated <- fundamentals
truncated[-exclind] <- NA
pipedata$actual.freq <- truncated

pipedata$actual.freq.full <- fundamentals

#semitone difference from reference frequencies
pipedata$out.of.tune <- 1200*log2(pipedata$actual.freq/pipedata$freq.eq.Hz)/100

#normalize to A = 440 Hz
pipedata$out.of.tune <- pipedata$out.of.tune - mean(pipedata$out.of.tune, na.rm = T)


#compare with equal temperament

#root mean square error
tuningrmse <- sqrt(mean(pipedata$out.of.tune^2, na.rm = T))


pipedata$rms <- rms/max(rms, na.rm = T)
pipedata$fundamp.norm <- fundamps/max(fundamps, na.rm = T)
pipedata$ampvar <- ampvar

pipedata$scentroid <- scentroid

pdf(width = 9, height = 5, 'temperament_vol.pdf', family = font)

function(){

print(
xyplot(freq.eq.Hz ~ key | stop, data = pipedata
, type = 'l', col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'frequency (Hz)'
, par.settings = list(strip.background = list(col = 'transparent'))
) +
xyplot(actual.freq ~ key | stop, data = pipedata
, col = 'red'
)
)

} #end of comment-out


#how many semitones out of tune
print(
barchart(out.of.tune ~ key | stop
, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'tuning error (semitone)'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)


#tuning error distribution
print(
densityplot(~ out.of.tune | stop, data = pipedata
, type = 'density'
, col = 'black'
, aspect = 1
, xlab = 'tuning error (semitone)', ylab = 'probability density'
, par.settings = list(strip.background = list(col = 'transparent'))
)
)



#stationary RMS volume

print(
barchart(rms ~ key | stop
, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'normalized RMS sound pressure'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)



#stationary fundamental volume

print(
barchart(fundamp.norm ~ key | stop, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'normalized fundamental sound pressure'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)




#time series amplitude variance

print(
barchart(ampvar ~ key | stop
, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'normalized sound pressure variance'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)



#spectral centroids

print(
barchart(scentroid ~ key | stop
, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'spectral centroid (Hz)'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)





#scatterplots of correlations between variables



panel = function (x, y, ...)
{
    panel.xyplot(x, y, ...)
    panel.text(x, y, labels = pipedata$pipe.id, pos = 3)
}

print(
xyplot(scentroid ~ actual.freq.full | stop, pipedata
, col = 'black'
, aspect = 1
, xlab = 'fundamental frequency (Hz)'
, ylab = 'spectral centroid (Hz)'
, par.settings = list(strip.background = list(col = 'transparent'))
#, panel = panel #label points
)
)




print(
xyplot(fundamp.norm ~ rms | stop, pipedata
, col = 'black'
, aspect = 1
, xlab = 'normalized RMS sound pressure'
, ylab = 'normalized fundamental sound pressure'
, par.settings = list(strip.background = list(col = 'transparent'))
)
)



print(
xyplot(ampvar ~ rms | stop, pipedata
, col = 'black'
, aspect = 1
, xlab = 'normalized RMS sound pressure'
, ylab = 'normalized sound pressure variance'
, par.settings = list(strip.background = list(col = 'transparent'))
)
)



print(
xyplot(ampvar ~ fundamp.norm | stop, pipedata
, col = 'black'
, aspect = 1
, xlab = 'normalized fundamental sound pressure'
, ylab = 'normalized sound pressure variance'
, par.settings = list(strip.background = list(col = 'transparent'))
)
)





dev.off()
embedFonts('temperament_vol.pdf')








#====================================
#attack transients
#====================================



ind1 <- samprate*attacktime[1]
ind2 <- samprate*attacktime[2]


print(noquote(paste('plotting attack transient spectrograms in',
folder, 'folder')))

pdf(height = 5, width = 5, 'attack_spectrograms.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

#power spectrogram
spectrograms <- list(NULL)
for(i in clipid)
    {
        print(noquote(paste('pipe', i)))
        
        spectrogram <- powspec(sounddata[[i]]@left[ind1:ind2], sr = samprate)
        spectrograms[[i]] <- spectrogram
        
        #frequency cutoff at specified harmonic
        afreqcutoff <- fundamentals[i]*ahcutoff
        cutoffrow <- min(1, afreqcutoff/(samprate/2))
        cutoffrow <- floor(nrow(spectrogram)*cutoffrow)
        
        image(seq(attacktime[1], attacktime[2], length.out = ncol(spectrogram))
        , seq(0, samprate/2, length.out = nrow(spectrogram))[1:cutoffrow]
        , log10(t(spectrogram))[,1:cutoffrow]
        , col = timPalette(10000)
        , mgp = mgp, tck = -0.01, axes = F
        , xlab = 'time (s)', ylab = 'frequency (Hz)'
        )

        axis(1, mgp = mgp, tck = -0.01)
        axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')
        
        mtext(paste('pipe', i, 'to harmonic', ahcutoff - 1))
        
    }


dev.off()
embedFonts('attack_spectrograms.pdf')


save(spectrograms, file = 'spectrograms.RData')






print(noquote(paste('plotting attack transient cepstrograms in',
folder, 'folder')))

pdf(height = 5, width = 5, 'attack_cepstrograms.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

#power spectrogram
cepstrograms <- list(NULL)
for(i in clipid)
    {
        print(noquote(paste('pipe', i)))
        
        
        subdata <- sounddata[[i]]
        subdata@left <- subdata@left[ind1:ind2]

        cepstrum <- melfcc(subdata
        , sr = samprate
        , maxfreq = 20000
        , numcep = numcep
        )

        cepstrograms[[i]] <- cepstrum
        
        image(cepstrum[,-1]
        , col = timPalette(10000)
        , xlab = 'time', ylab = 'melfcc'
        , mgp = mgp, tck = -0.01, axes = F
        )

        axis(1, mgp = mgp, tck = -0.01)
        axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')
        
        mtext(paste('pipe', i, 'attack'))
        
    }


dev.off()
embedFonts('attack_cepstrograms.pdf')

save(cepstrograms, file = 'cepstrograms.RData')









save(pipedata, file = 'pipedata.RData')