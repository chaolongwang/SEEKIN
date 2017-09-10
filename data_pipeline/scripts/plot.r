args<-commandArgs()

input=args[1+5]
output=args[2+5]
kin <- read.table(input, header=T);

n=dim(kin)[1] 
p.cols <- rep("white", n);
p.cols[which(kin$relationship=="MZ_Twin")] <- "purple";
p.cols[which(kin$relationship=="PO")] <- "orange";
p.cols[which(kin$relationship=="FS")] <- "red";
p.cols[which(kin$relationship=="2nd")] <- "blue";
p.cols[which(kin$relationship=="3nd")] <- "green";
p.cols[which(kin$relationship=="Unrelated")] <- "black";
p.cols[which(kin$relationship=="Ambiguous")] <- "white";

removeList<-c("ME-G1Y4XFM","ME-YIPXCPE","ME-JNIXHIA","ME-2VQVD3C","ME-E5JBYEP","ME-2PUIHX4","ME-E5JBYEP","ME-KS447JZ","ME-HHC8NU7","ME-YIPXCPE")
pos1=which(kin$Ind1%in%removeList);
pos2=which(kin$Ind2%in%removeList);
pos_remove=union(pos1,pos2);


mainVector=c("Pairs of Merck dataset","Pairs of Malays","Pairs of Chinese and Malay");
leg.txt <- c("Monozygotic twin", "Parent-offspring","Full sibling", "2nd degree", "3rd degree", "Unrelated")
leg.cols <- c("purple","orange","red","blue","green","black");
Type<-c("PO","FS","2nd","3nd","Unrelated")
Type<-c("Unrelated","PO","2nd","3nd")

PopList=c(11,22,12)

png(paste(output,".png",sep=''),  width=6, height=6.5, units="in", res=300)
methodList<-c(4,4,4)  #...Change lcMLkin
par(mfrow=c(1,1))
RMSE <- matrix(1:30,nrow=3,ncol=10)
RMSE1<- matrix(1:30,nrow=3,ncol=10)
mymin1=min(kin[,4])
mymax1=max(kin[,4])
mymin2=min(kin$Kin_true) 	# Change the threshold if the figure is not good
mymax2=max(kin$Kin_true)
mymin=min(mymin1,mymin2)
mymax=max(mymax1,mymax2)

True=kin$Kin_true;

for(i in 1:1){
	 Pos=which(kin$relationship!="Ambiguous" )
	 #Pos=which(kin$relationship!="Ambiguous" & (kin$class==PopList[i]))
     #if(i==3){
	 #	Pos=which(kin$relationship!="Ambiguous" & (kin$class==12 | kin$class==21))
	 #}
	 Pos=setdiff(Pos,pos_remove)


 	for(j in 1:4){
        if(j==2){  RelationPos=which(kin$relationship=="PO" | kin$relationship=="FS")}
        else{ RelationPos=which(kin$relationship==Type[j]) }
        RelationPos=intersect(RelationPos, Pos)
        Num=length(RelationPos)
		print(Num)
        x<-True[RelationPos]-kin[RelationPos,methodList[i]];
        RMSE[i,2*j-1]<-round(sqrt(sum(x^2))/sqrt(Num),4);
        RMSE[i,2*j]<-round(mean(x),4);
    }

    # All 
    j=5;
    x<-True[Pos]-kin[Pos,methodList[i]];
    RMSE[i,2*j-1]<-round(sqrt(sum(x^2))/sqrt(length(Pos)),4);
    RMSE[i,2*j]<-round(mean(x),4)
 

	
    plot(kin$Kin_true[Pos],kin[Pos,methodList[i]],col=p.cols[Pos],xlab="Array",ylab="Seq",main=mainVector[i], asp=1, xlim=c(mymin, mymax), ylim=c(mymin, mymax))
    abline(a=0, b=1, col="black")
    if(i==1){
      legend("topleft", leg.txt, fill=leg.cols, ncol=1, cex=0.6);
    }
}
dev.off()



