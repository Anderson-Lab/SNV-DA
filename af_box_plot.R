library(ggplot2)

exonic = read.csv("nonsyn_for_boxplot.csv", header=F)

df = data.frame(snp = factor(exonic[,1], levels=unique(exonic[,1])), af = exonic[,2], Outcome=exonic[,3])

pdf("nonsyn_af_boxplot.pdf", 7,3)
plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=8)
plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
plot4 = plot4 + theme_bw(base_size=10.5) 
plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
plot4 = plot4 + ylab("Allele Fraction")
plot4 = plot4 + xlab(NULL) + ggtitle("Top Nonsynonymous Exonic SNVs")
plot4 = plot4 + theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(), legend.position="bottom")
#plot4 = plot4 + theme(legend.position="bottomright",plot.background = element_blank(),
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#panel.border = element_blank(),axis.line = element_line(color='black'))
plot4
dev.off()

exonic = read.csv("intronic_for_boxplot.csv", header=F)

df = data.frame(snp = factor(exonic[,1], levels=unique(exonic[,1])), af = exonic[,2], Outcome=exonic[,3])

pdf("intronic_af_boxplot.pdf", 7,3)
plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=8)
plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
plot4 = plot4 + theme_bw(base_size=10.5) 
plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
plot4 = plot4 + ylab("Allele Fraction")
plot4 = plot4 + xlab(NULL) + ggtitle("Top Intronic SNVs")
plot4 = plot4 + theme(axis.title.x=element_blank(),
											axis.text.x=element_blank(),
											axis.ticks.x=element_blank(), legend.position="bottom")
#plot4 = plot4 + theme(legend.position="bottomright",plot.background = element_blank(),
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#panel.border = element_blank(),axis.line = element_line(color='black'))
plot4
dev.off()





exonic = read.csv("pulled_top_UTR3.csv", header=F)

df = data.frame(snp = factor(exonic[,1], levels=unique(exonic[,1])), af = exonic[,2], Outcome=exonic[,3])

pdf("UTR3_plot3.pdf", 7,3)
plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=8)
plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
plot4 = plot4 + theme_bw(base_size=10.5) 
plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
plot4 = plot4 + ylab("Allele Fraction")
plot4 = plot4 + xlab(NULL) + ggtitle("Top 3'UTR SNV Features")
plot4 = plot4 + theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(), legend.position="bottom")
#plot4 = plot4 + theme(legend.position="bottomright",plot.background = element_blank(),
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#panel.border = element_blank(),axis.line = element_line(color='black'))
plot4
dev.off()


exonic = read.csv("pulled_top_coding.csv", header=F)

df = data.frame(snp = factor(exonic[,1], levels=unique(exonic[,1])), af = exonic[,2], Outcome=exonic[,3])

pdf("coding_plot3.pdf", 7,3)
plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=8)
plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
plot4 = plot4 + theme_bw(base_size=10.5) 
plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
plot4 = plot4 + ylab("Allele Fraction")
plot4 = plot4 + xlab(NULL) + ggtitle("Top Exonic + UTRs SNV Features")
plot4 = plot4 + theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(), legend.position="bottom")
#plot4 = plot4 + theme(legend.position="bottomright",plot.background = element_blank(),
# panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#panel.border = element_blank(),axis.line = element_line(color='black'))
plot4
dev.off()


exonic = read.csv("pulled_top_ab.csv", header=F)

df = data.frame(snp = factor(exonic[,1], levels=unique(exonic[,1])), af = exonic[,2], Outcome=exonic[,3])
df
pdf("ab_plot3.pdf", 7,3)
plot4 = ggplot(df, aes(x=Outcome,y=af)) + geom_boxplot( frame.plot=FALSE,axes=FALSE,outlier.size=0) + facet_grid(. ~ snp) + facet_wrap( ~ snp, ncol=8)
plot4 = plot4 +  geom_point(data=subset(df,af>0),size=1.25, aes(colour=Outcome),position = position_jitter(width = 0.1)) + scale_color_brewer(palette="Set1")
plot4 = plot4 + theme_bw(base_size=10.5) 
plot4 = plot4 + scale_y_continuous(breaks=c(0,.5,1))
plot4 = plot4 + ylab("Allele Fraction")
plot4 = plot4 + xlab(NULL) + ggtitle("Top IG SNV Features")
plot4 = plot4 + theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(), legend.position="bottom")
#plot4 = plot4 + theme(legend.position="bottomright",plot.background = element_blank(),
                     # panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      #panel.border = element_blank(),axis.line = element_line(color='black'))
plot4
dev.off()
#pdf("boxplot_combined2.pdf",7,10)
#multiplot(plot1, plot2, plot3, plot4, cols=1)


