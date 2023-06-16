require (ggplot2)
require(sjPlot)


dados=read.table(pipe("pbpaste"), sep="\t", header=T, row.names=1);dados
names(dados)


dados1=as.data.frame(scale(dados))
mod=lm(IVI ~  Habitat.Diversity+ Distance.Nearest.Patch+ Area, data=dados1)

plot1=plot_model(mod, colors="bw", show.values=TRUE, title= "Juncais", axis.title ="")
plot1 <- plot1 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot2=plot_model(mod, colors="bw", show.values=TRUE, title= "Sangradouros", axis.title ="")

plot2 <- plot2 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

grid.arrange(plot1,plot2, nrow=1)

