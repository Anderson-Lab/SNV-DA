library("ggplot2")

input = read.table("nonsyn_parallel_iteration_optimals.txt", header=F)
input
plot = plot(hist(unlist(input), breaks=50))
plot = plot + plot(density(unlist(input)))

hist(unlist(input), breaks=35, freq=F)
lines(density(unlist(input)), lwd=2)


plot = ggplot(unlist(input))
plot = plot + geom-histrogram()
plot = plot + geom_histogram(aes(y=..density..)) + geom_density()
plot


plot = qplot(unlist(input), geom='histogram', binwidth=10, frequency=F
plot = plot + geom_density()
plot

pdf("nonsyn_density.pdf", 7,3)
d = density(unlist(input))
optM = round(d$x[which.max(d$y)])
dataset <- data.frame(X = unlist(input))
plot = ggplot(dataset, aes(x = X)) + geom_histogram(aes(y = ..density..), colour='black', fill='blue',binwidth=35) + geom_density(color='red', size=1)
plot = plot + geom_vline(xintercept=c(optM), linetype="dashed", size=1)
plot = plot + xlab("Value of K Tested")
plot = plot + ylab("Density")
plot = plot + scale_y_continuous(expand = c(0,0), limits=c(0,.0052)) + scale_x_continuous(limits = c(0,max(unlist(input))),expand = c(0,0))
plot = plot + theme_bw()
plot = plot + theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold"))
#plot = plot + theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(color='black'))
plot
dev.off()