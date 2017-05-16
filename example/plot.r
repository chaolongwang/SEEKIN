args<-commandArgs()

inputEst=args[1+5]
inputTrue=args[2+5]
output=args[3+5]
kinEst <- read.table(inputEst, header=T);
kin <- read.table(inputTrue, header=T);



n=dim(kin)[1] 
p.cols <- rep("white", n);
p.cols[which(kin$relationship=="MZ_Twin")] <- "purple";
p.cols[which(kin$relationship=="PO")] <- "orange";
p.cols[which(kin$relationship=="FS")] <- "red";
p.cols[which(kin$relationship=="2nd")] <- "blue";
p.cols[which(kin$relationship=="3nd")] <- "green";
p.cols[which(kin$relationship=="Unrelated")] <- "black";
p.cols[which(kin$relationship=="Ambiguous")] <- "white";


leg.txt <- c("Monozygotic twin", "Parent-offspring","Full sibling", "2nd degree", "3rd degree", "Unrelated")
leg.cols <- c("purple","orange","red","blue","green","black");
Type<-c("PO","FS","2nd","3nd","Unrelated")
Type<-c("Unrelated","PO","2nd","3nd")


png(paste(output,".png",sep=''), width=5, height=5, units="in", res=300)
par(mfrow=c(1,1))
par(pty="s")
mymin=min(kinEst$kinship)
mymax=max(kinEst$kinship)
plot(kin$kin,kinEst$kinship,col=p.cols,xlab="Array",ylab="Seq",asp=1, xlim=c(mymin, mymax), ylim=c(mymin, mymax))
abline(a=0, b=1, col="black")
legend("topleft", leg.txt, fill=leg.cols, ncol=1, cex=0.8);

dev.off()



