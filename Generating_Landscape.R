## This code is in base R and was used for my preprint "A Waddington Epigenetic Landscape for the C. elegans embryo" (bioRxiv 2020).
## The code converts the C. elegans early embryonic lineage to a matrix and uses the number of detected genes per cell (scRNAseq dataset Tintori et al. 2016)
## in order to define the elevation of each cell region in the matrix. The result is a downward slope communicating the decrease in the number of
## expressed genes as cells differentiate. The landscape can then be used to plot gene expression during early embryogenesis on a "Developmental Spacetime"
## that encompasses all cells from generated during the first four rounds of division.
## Note that in the preprint I use a different metric than "number of detected genes" to define elevation. There, I used "Epigenetic tension" (the formula is in the preprint)
## Currently, I prefer to use Shannon´s Entropy to define elevation. 

library(RColorBrewer)

#Assign Tintori dataset samples to C. elegans embryonic cells
P0 <- c('st451','st449','st441','st413','st411')
AB <- c('st409','st361','st311','st301','st265')
P1 <- c('st410','st362','st312','st266','st302')
EMS <- c('st406','st364','st335','st305','st193')
P2 <- c('st408','st366','st336','st308','st196')
MS <- c('st485','st399','st374','st342','st229')
MSx1 <- c('st513','st476','st529','st323','st257')
MSx2 <- c('st517','st474','st528','st354','st353','st324','st258')
E <- c('st486','st398','st375','st343','st231')
Ea <- c('st514','st475','st356','st260')
Ep <- c('st516','st473','st527','st355','st325','st261')
P3 <- c('st488','st400','st376','st344','st232')
C <- c('st487','st397','st373','st341','st230')
Cx1 <- c('st515','st531','st357','st322','st262')
Cx2 <- c('st518','st477','st530','st358','st321','st259')
D <- c('st520','st479','st533','st359','st327','st264')
P4 <- c('st519','st478','st532','st360','st328','st263')
ABa <- c('st407','st363','st334','st306','st194')
ABp <- c('st405','st365','st333','st307','st195')
ABar <- c('st482','st393','st372','st338','st227')
ABal <- c('st484','st396','st370','st226','st337')
ABpr <- c('st481','st395','st369','st339','st228')
ABpl <- c('st483','st394','st371','st340','st225')
ABplx <-c('st510','st507','st471','st468','st542','st541','st349','st347','st317','st252','st249')
ABprx <- c('st509','st505','st472','st465','st540','st537','st352','st350','st318','st313','st251','st250') 
ABalx <- c('st469','st467','st539','st536','st351','st345','st319','st314','st256','st255')
ABarx <- c('st508','st506','st470','st466','st538','st535','st348','st346','st320','st316','st254','st253')


## Number of detected genes in R
rpkm <- read.csv('TableS2_RPKMs.csv',header=TRUE)
rownames(rpkm) <- rpkm$transcript
rpkm <- rpkm[,2:ncol(rpkm)]

ndg <- list()
for (i in colnames(rpkm)){
  ndg[[i]] <- length(rpkm[rpkm[,i]>0,i]) 
  }

ndg <- as.matrix(ndg)
hl = max(as.numeric(ndg))+100 # high limit
ll = min(as.numeric(ndg))-100 # low limit

# Convert embryonic lineage to landscape with elevation determined by number of detected genes

# Generate a lineage matrix called ndgheat (i.e. Number of Detected Genes Heatmap)
ndgheat <- matrix(0L,16,32)

# Populate regions of the matrix with data from it's corresponding cell
ndgheat[1,1:32]                                   <- mean(unlist(ndg[P0,]))#P0 region
ndgheat[2:5,  ((32/2)*1)+1  :  as.integer(32/2) ]	<- mean(unlist(ndg[P1,]))#P1 region
ndgheat[6:10, ((32/4)*3)+1  :  as.integer(32/4) ]	<- mean(unlist(ndg[P2,]))#P2 region
ndgheat[11:15,((32/8)*7)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[P3,]))#P3 region
ndgheat[16,   ((32/16)*15)+1:  as.integer(32/16)]	<- mean(unlist(ndg[P4,]))#P4 region
ndgheat[2:4,  1             :  as.integer(32/2) ]	<- mean(unlist(ndg[AB,]))#AB region
ndgheat[6:9,  ((32/4)*2)+1  :  as.integer(32/4) ]	<- mean(unlist(ndg[EMS,]))#EMS region
ndgheat[10:15,((32/8)*5)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[E,]))#E region
ndgheat[16,   ((32/16)*10)+1:  as.integer(32/16)]	<- mean(unlist(ndg[Ea,]))#Ea region
ndgheat[16,   ((32/16)*11)+1:  as.integer(32/16)]	<- mean(unlist(ndg[Ep,]))#Ep region
ndgheat[10:15,((32/8)*4)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[MS,]))#MS region
ndgheat[16,   ((32/16)*8)+1 :  as.integer(32/16)]	<- mean(unlist(ndg[MSx1,]))#MSx1 region
ndgheat[16,   ((32/16)*9)+1 :  as.integer(32/16)]	<- mean(unlist(ndg[MSx2,]))#MSx2 region
ndgheat[11:15,((32/8)*6)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[C,]))#C region
ndgheat[16,   ((32/16)*12)+1:  as.integer(32/16)]	<- mean(unlist(ndg[Cx1,]))#Cx1 region
ndgheat[16,   ((32/16)*13)+1:  as.integer(32/16)]	<- mean(unlist(ndg[Cx2,]))#Cx2 region
ndgheat[16,   ((32/16)*14)+1:  as.integer(32/16)]	<- mean(unlist(ndg[D,]))#D region
ndgheat[5:8,  1             :  as.integer(32/4) ]	<- mean(unlist(ndg[ABa,]))#ABa region
ndgheat[5:8,  ((32/4)*1)+1  :  as.integer(32/4) ]	<- mean(unlist(ndg[ABp,]))#ABp region
ndgheat[9:11, ((32/8)*1)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[ABar,]))#ABar region
ndgheat[9:11, 1             :  as.integer(32/8) ]	<- mean(unlist(ndg[ABal,]))#ABal region
ndgheat[9:11, ((32/8)*3)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[ABpr,]))#ABpr region
ndgheat[9:11, ((32/8)*2)+1  :  as.integer(32/8) ]	<- mean(unlist(ndg[ABpl,]))#ABpl region
ndgheat[12:16,1             :  as.integer(32/8) ]	<- mean(unlist(ndg[ABalx,]))#ABalx region
ndgheat[12:16,((32/16)*2)+1 :  as.integer(32/8) ]	<- mean(unlist(ndg[ABarx,]))#ABarx region
ndgheat[12:16,((32/16)*4)+1 :  as.integer(32/8) ]	<- mean(unlist(ndg[ABplx,]))#ABplx region
ndgheat[12:16,((32/16)*6)+1 :  as.integer(32/8) ]	<- mean(unlist(ndg[ABprx,]))#ABprx region

# Represent number of detected genes by elevation of terrain (z axis) and with terrain colors
col.pal.ndg<- terrain.colors
colors.ndg<-col.pal.ndg(100)
z<-ndgheat
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
z.facet.range<-cut(z.facet.center, 100)

# Plot landscape
par(mfrow=c(2,2))
persp(as.matrix(ndgheat),expand = 0.3,theta = 120,phi=30,ticktype = 'detailed',cex.axis=0.5,col=colors.ndg[z.facet.range],
      xlab='Distance from Zygote',ylab='',zlab='Detected Genes',cex.lab=0.5)
persp(as.matrix(ndgheat),expand = 0.3,theta=90,phi=30,ticktype = 'detailed',cex.axis=0.5,col=colors.ndg[z.facet.range],
      xlab='',ylab='',zlab='Detected Genes',cex.lab=0.5)
persp(as.matrix(ndgheat),expand = 0.3,theta=0,phi=0,ticktype = 'detailed',cex.axis=0.5,col=colors.ndg[z.facet.range],
      xlab='Distance from Zygote',ylab='',zlab='Detected Genes',cex.lab=0.5)
persp(as.matrix(ndgheat),expand = 0.3,theta=180,phi=0,ticktype = 'detailed',cex.axis=0.5,col=colors.ndg[z.facet.range],
      xlab='Distance from Zygote',ylab='',zlab='Detected Genes',cex.lab=0.5)


# Use landscape to plot gene expression. The first three genes represent maternal regulators
# The second three genes represent mesoendoderm master regulators
# The final trio represent RNAi pathway genes
col.pal<-colorRampPalette(rev(brewer.pal(n = 3, name = "RdBu")))
colors<-col.pal(100)
par(mfrow=c(3,3))
for (i in c('pie-1','pos-1','F32D1.6','skn-1','med-1','end-1','dcr-1','rde-1','alg-1')){
  irpkm <- matrix(0L,16,32)
  irpkm[1,1:32]			<- mean(as.numeric(rpkm[i,P0]))	#P0 region
  irpkm[2:5,  ((32/2)*1)+1  :  as.integer(32/2) ]	<- mean(as.numeric(rpkm[i,P1]))	#P1 region
  irpkm[6:10, ((32/4)*3)+1  :  as.integer(32/4) ]	<- mean(as.numeric(rpkm[i,P2]))	#P2 region
  irpkm[11:15,((32/8)*7)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,P3]))	#P3 region
  irpkm[16,   ((32/16)*15)+1:  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,P4]))	#P4 region
  irpkm[2:4,  1             :  as.integer(32/2) ]	<- mean(as.numeric(rpkm[i,AB]))	#AB region
  irpkm[6:9,  ((32/4)*2)+1  :  as.integer(32/4) ]	<- mean(as.numeric(rpkm[i,EMS]))	#EMS region
  irpkm[10:15,((32/8)*5)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,E]))	#E region
  irpkm[16,   ((32/16)*10)+1:  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,Ea]))	#Ea region
  irpkm[16,   ((32/16)*11)+1:  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,Ep]))	#Ep region
  irpkm[10:15,((32/8)*4)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,MS]))	#MS region
  irpkm[16,   ((32/16)*8)+1 :  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,MSx1]))	#MSx1 region
  irpkm[16,   ((32/16)*9)+1 :  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,MSx2]))	#MSx2 region
  irpkm[11:15,((32/8)*6)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,C]))	#C region
  irpkm[16,   ((32/16)*12)+1:  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,Cx1]))	#Cx1 region
  irpkm[16,   ((32/16)*13)+1:  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,Cx2]))	#Cx2 region
  irpkm[16,   ((32/16)*14)+1:  as.integer(32/16)]	<- mean(as.numeric(rpkm[i,D]))	#D region
  irpkm[5:8,  1             :  as.integer(32/4) ]	<- mean(as.numeric(rpkm[i,ABa]))	#ABa region
  irpkm[5:8,  ((32/4)*1)+1  :  as.integer(32/4) ]	<- mean(as.numeric(rpkm[i,ABp]))	#ABp region
  irpkm[9:11, ((32/8)*1)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,ABar]))	#ABar region
  irpkm[9:11, 1             :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,ABal]))	#ABal region
  irpkm[9:11, ((32/8)*3)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,ABpr]))	#ABpr region
  irpkm[9:11, ((32/8)*2)+1  :  as.integer(32/8) ]	<- mean(as.numeric(rpkm[i,ABpl]))	#ABpl region
  irpkm[12:16,1             :  as.integer(32/8)]	<- mean(as.numeric(rpkm[i,ABalx]))	#ABalx region
  irpkm[12:16,((32/16)*2)+1 :  as.integer(32/8)]	<- mean(as.numeric(rpkm[i,ABarx]))	#ABarx region
  irpkm[12:16,((32/16)*4)+1 :  as.integer(32/8)]	<- mean(as.numeric(rpkm[i,ABplx]))	#ABplx region
  irpkm[12:16,((32/16)*6)+1 :  as.integer(32/8)]	<- mean(as.numeric(rpkm[i,ABprx]))	#ABprx region
  z<-irpkm
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)
  persp(as.matrix(ndgheat),expand = 0.3, theta = 90,phi=30,ticktype = 'detailed',cex.axis=0.5, col=colors[z.facet.range],xlab='',ylab='',zlab='',border='black',lwd=0.05, box=FALSE, main=i)
  legend(0.24,0.1,legend="",fill=colors[1],bty='n',border=NA)
  legend(0.24,0.13,legend="",fill=colors[50],bty='n',border=NA)
  legend(0.24,0.16,legend="",fill=colors[100],bty='n',border=NA)
  text(0.37,0.09,round(log(max(irpkm))),cex=0.55)
  text(0.37,0.06,round(log(max(irpkm))/2),cex=0.55)
  text(0.37,0.03,'0',cex=0.55)
  text(0.34,0.13,'RPKM',cex=0.55)
  text(0.34,0.16,'ln',cex=0.55)
}

