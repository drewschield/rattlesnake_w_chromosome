setwd('./')

bsqual <- read.table('./resources/fastqc_per_base_sequence_quality_plot.tsv',header=T,sep='\t')
psqual <- read.table('./resources/fastqc_per_sequence_quality_scores_plot.tsv',header=T,sep='\t')
nuccomp <- read.table('./resources/multiqc_nuc_comp.txt',header=T,sep='\t')
psgc <- read.table('./resources/fastqc_per_sequence_gc_content_plot.tsv',header=T,sep='\t')

par(mfrow=c(2,2))
plot(bsqual$Position..bp.,((bsqual$PitViper_042019_S2_L004_R1_001+bsqual$PitViper_042019_S2_L004_R2_001)/2),type='l',lwd=2,ylim=c(0,38),xlab='Position (bp)',ylab='Phred Score',main='Mean Quality Scores Per Base')
abline(h=28,lty=2)
abline(h=20,lty=2)

plot(psqual$Mean.Sequence.Quality..Phred.Score.,((psqual$PitViper_042019_S2_L004_R1_001+psqual$PitViper_042019_S2_L004_R2_001)/2),type='l',lwd=2,xlab='Phred Score',ylab='Read Count',main='Per Sequence Quality Scores')
abline(v=28,lty=2)
abline(v=20,lty=2)

plot(nuccomp$Position..bp.,nuccomp$X..T,type='l',lwd=2,col='red',ylim=c(0,100),xlab='Position (bp',ylab='% Reads',main='Per Base Sequence Content')
lines(nuccomp$Position..bp.,nuccomp$X..C,lwd=2,col='blue')
lines(nuccomp$Position..bp.,nuccomp$X..A,lwd=2,col='green')
lines(nuccomp$Position..bp.,nuccomp$X..G,lwd=2,col='black')
legend(x=0,y=100,legend=c('T','C','A','G'),col=c('red','blue','green','black'),lwd=2)

plot(psgc$X..GC,((psgc$PitViper_042019_S2_L004_R1_001+psgc$PitViper_042019_S2_L004_R2_001)/2),type='l',lwd=2,xlab='% GC',ylab='% Reads',main='Per Sequence GC Content')
