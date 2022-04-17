############################ Import Libraries ############################
############################ You might need to install these with "install.packages" or with Bioconductor https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(ggtree)
library(ggmsa)
library(ggplot2)
library(treeio)
library(ape)
library(phytools)
library(Biostrings)
library(grid)
library(gtable)

############################ Set Working Directory. You'll have to change this on your computer ############################ 
setwd('/home/sdenecke/Dropbox/Shane/target_discovery/Lep_Target_Characterization/GGPPS/OrthoDB_analysis')

name.short=function(x){
  x=gsub('(^.+)__.+$','\\1',x)
  x=gsub('_GGPPS_LOC118281562','',x)
  if(x=='Sc_2DH4_A'){
    return(x)
  }else{
    a=strsplit(x,'_') %>% unlist()
    b=substr(a[1],1,1)
    return(paste0(b,'. ',a[2]))
  }
}

############################  Read in Tree and make Yeast the outgroup ############################ 
tree=read.tree('./outputs/GGPPS_tree.nwk') %>% midpoint.root()
#tree=ape::root(tree,outgroup='Sc_2DH4_A')
tree$tip.label=gsub('_GGPPS_LOC118281562','',tree$tip.label)
tree$tip.label=sapply(tree$tip.label,name.short)

############################Import The Multiple Sequence Alignment ############################ 
strings=readAAStringSet('./outputs/GGPPS_named.faa.aln.trimm')

data1 = tidy_msa(strings, 33, 38) ###### Subset out AA 90-110 on MSA
data1$name=sapply(data1$name,name.short)

data2 = tidy_msa(strings, 54, 59) ###### Subset out AA 90-110 on MSA
data2$name=sapply(data2$name,name.short)


data3 = tidy_msa(strings, 98, 103) ###### Subset out AA 90-110 on MSA
data3$name=sapply(data3$name,name.short)


data4 = tidy_msa(strings, 161, 166) ###### Subset out AA 90-110 on MSA
data4$name=sapply(data4$name,name.short)


data5 = tidy_msa(strings, 185, 190) ###### Subset out AA 90-110 on MSA
data5$name=sapply(data5$name,name.short)



############################ Make Tree ############################
gp=ggtree(tree,size=1)
gp=gp+geom_tippoint(size=4)
gp=gp+xlim_tree(2)
gp=gp+labs(x='Substitutions / Site')
gp=gp+geom_tiplab(fontface='bold',hjust=-.3,align=T,size=8)
gp=gp+theme(legend.text=element_text(size=30,face='bold'),legend.title = element_text(size=36,face='bold'),strip.text = element_text(size= 40,face='bold'))
gp=gp+theme_tree2()
############################ Print only the tree ############################
print(gp)

############################ Add the MSA ############################
gp2=facet_plot(gp,geom=geom_msa,posHighligthed = c(36),panel='Region 1',data=data1,char_width=.4)
gp3=facet_plot(gp2,geom=geom_msa,posHighligthed = c(57),panel='Region 2',data=data2,char_width=.4)
gp4=facet_plot(gp3,geom=geom_msa,posHighligthed = c(100),panel='Region 3',data=data3,char_width=.4)
gp5=facet_plot(gp4,geom=geom_msa,posHighligthed = c(163),panel='Region 4',data=data4,char_width=.4)
gp6=facet_plot(gp5,geom=geom_msa,posHighligthed = c(188),panel='Region 5',data=data5,char_width=.4)

gp7=gp4+theme(legend.position = 'none',strip.text=element_text(size=36,face='bold'),axis.text=element_text(size=14))
print(gp7)


############## Adjust panel sizes
gt = ggplot_gtable(ggplot_build(gp7))
#gtable_show_layout(gt) # will show you the layout - very handy function
#gt # see plot layout in table format
gt$layout$l[grep('panel-1', gt$layout$name)] # you want to find the column specific to panel-2
#gt$widths[5] = 0.8*gt$widths[5] # in this case it was colmun 7 - reduce the width by a half
gt$widths[7] = 0.25*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
gt$widths[9] = 0.25*gt$widths[9] # in this case it was colmun 7 - reduce the width by a half
gt$widths[11] = 0.25*gt$widths[11] # in this case it was colmun 7 - reduce the width by a half
gt$widths[13] = 0.25*gt$widths[13] # in this case it was colmun 7 - reduce the width by a half
gt$widths[15] = 0.25*gt$widths[15] # in this case it was colmun 7 - reduce the width by a half

grid.draw(gt) # plot with grid draw

ggsave(filename='tree_msa_scaled.pdf',plot=gt,device='pdf',height=10,width=25)


# 
# gp7=facet_plot(gp,geom=geom_msa,panel='Multiple Sequence Alignment',data=data6,font = NULL, color = "Chemistry_AA")
# 
# 
# gp8=facet_plot(gp,geom=geom_msa,panel='Region 1',data=data7,posHighligthed=c(36,57),char_width=.7)
# gp9=facet_plot(gp8,geom=geom_msa,panel='Region 2',data=data8,posHighligthed=c(100),char_width=.7)
# gp10=facet_plot(gp9,geom=geom_msa,panel='Region 3',data=data9,posHighligthed=c(163,188),char_width=.7)
# 
# gp10=gp10+theme(legend.position = 'none',strip.text=element_text(size=20))
# 
# ############################ Save plot ############################
# ggsave(filename='Plot_output.pdf',plot=gp10,device='pdf',height=10,width=20)


