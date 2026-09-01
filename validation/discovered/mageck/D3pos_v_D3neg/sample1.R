pdf(file='sample1.pdf',width=4.5,height=4.5);
gstable=read.table('sample1.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("SMAD5","SMAD1","DBR1","ACVR2A","DROSHA","CDH5","BMPR1A","NAA25","TBX3","ASCC3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='D3PR1,D3PR2_vs_D3NR1,D3NR2 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(16460.35888878052,13182.524887236528,1359.9387325554028,2420.600617321856),c(6414.076800710754,6606.912276433313,208.95345442592617,682.4434762231401),c(19541.675897817615,15711.537870148468,980.3399570149702,2151.8245363152682),c(11053.196140355256,12430.43863595331,1911.9241079972244,1030.936711516251),c(16882.097010566384,15524.634152528737,1370.386405276699,1864.196425578759),c(14838.552940774716,13068.057538646453,969.8922842936739,1297.5582748394195))
targetgene="SMAD5"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(7864.444487827993,9404.208107603237,212.4360119996916,1247.4657386999152),c(5428.878261471147,8132.547406860388,799.2469631791676,652.8188580761214),c(8197.034713951803,12736.281082967413,644.2731511466056,1120.3491953781622),c(4256.24055699338,5488.172799197515,142.78486052438288,628.5805341376515),c(6265.496974607332,7604.030195792157,860.1917207200627,659.28241112638),c(9363.957809733287,11792.819733260165,600.7411814745377,1019.0868642574436))
targetgene="SMAD1"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(14859.125532081343,11144.4695166367,1481.828247637193,884.4295090437222),c(15709.459306088615,12003.868907223114,471.8865512452166,1811.9493717558353),c(5928.335061526493,4525.037373951349,557.2092118024698,831.6438257999434),c(9705.691409771154,11359.990071403949,994.2701873100319,985.1532107435859),c(8935.36215751188,6459.356709891421,583.3283936057105,477.76429629828357),c(6491.795478980235,4753.972071131496,99.25289085231492,812.7917960700224))
targetgene="DBR1"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(14782.549775551119,15208.060391584315,799.2469631791676,1657.362727970483),c(7635.860139976576,7644.272623030855,296.01739377006203,946.9105218628889),c(15419.157184317315,13924.774100750286,585.0696723925932,1311.5626397816466),c(3355.6182264587974,3555.6420155791625,5491.993293828093,2913.4465374040765),c(9974.278018496569,13065.374710163875,306.46506649135836,986.2304695852956),c(11399.501427350153,8997.31245441196,1290.287581080094,611.3443926702952))
targetgene="ACVR2A"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(6463.2224354988075,9575.014854327488,677.3574480973773,1536.7097376989886),c(4059.658017841161,4203.097956041766,532.8313087861117,449.7555664138295),c(3569.344591699872,6737.476595918865,391.78772704861154,748.694894988291),c(5164.863339702761,7137.218039823263,933.3254297691368,736.8450477294836),c(12823.581914464477,14067.858286487879,2549.2321439962993,2402.28721701279),c(2864.1618785782516,4301.468333736361,538.0551451467599,279.01004000283075))
targetgene="DROSHA"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(7758.1527660770835,4742.346481040317,1694.2642596368846,1618.5814096689314),c(11802.952801307903,13669.01111874434,2439.531580422688,2207.841996084176),c(13475.047305841015,13348.86025315648,3783.798803896146,1894.3596731466328),c(7629.002609541033,7017.38503426803,1206.7061993097236,2141.0519478981705),c(12330.982644844675,11915.335567297978,2352.467641078552,2052.716722877969),c(5798.041983251185,6390.497445505205,1946.7496837348788,1353.5757346083276))
targetgene="CDH5"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(6783.240522490792,4974.858282863904,240.29647258981507,189.59755614091975),c(2134.977808932232,3535.967940040244,2063.415362456021,2655.981674235441),c(7350.129705162305,10226.047899433532,212.4360119996916,747.6176361465813),c(9562.826192364018,11951.106613732376,947.2556600641985,1024.4731584659926),c(8365.044209622594,9212.833009179209,182.83427262268538,624.2714987708124),c(14513.963166825704,17014.498236521416,1864.9095807513909,1616.9655214063666))
targetgene="BMPR1A"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5442.593322342233,5220.78422710039,1191.034690227779,791.7852486566819),c(13187.031027548232,10274.33881211997,2392.5170531768545,2136.7429125313315),c(2556.715930718096,2800.8729358133646,818.4010298348775,619.9624634039734),c(5616.317426709309,5474.758656784616,959.4446115723777,1046.5569647210427),c(3380.7625047224533,3996.520162883118,1058.6975024246926,410.43561869142286),c(6420.934331146296,5192.167389952871,1361.6800113422855,1436.52466541998))
targetgene="NAA25"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3945.365843915453,5354.031375068523,524.1249148516981,1022.8572702034278),c(2585.288974199523,3503.7739982492853,1191.034690227779,1205.9912732940888),c(5643.747548451479,4792.425946048474,940.2905449166677,1123.5809719032914),c(8991.365322735477,9863.866054285252,1323.3718780308657,1756.470541407782),c(6062.056905019572,7284.7736063651555,1932.819453439817,1564.7184675834428),c(5254.011235364814,8296.19994429776,583.3283936057105,959.8376279634062))
targetgene="TBX3"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5395.733531032692,5768.081237546679,1648.9910111779338,1798.4836362344631),c(1899.535930645273,2534.3786398770994,240.29647258981507,568.2540390019043),c(3159.0356873065793,3356.2184317073934,1379.0927992111126,948.5264101254536),c(3645.920348230097,5177.8589713791125,853.2266055725318,1099.3426479648215),c(2979.596974243217,2420.8055674478856,806.2120783266985,520.8546499666744),c(5723.7520701994745,5050.871756537,243.77903016358053,1098.8040185439668))
targetgene="ASCC3"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("NON-TARGETING","CBLL1","SUZ12","SAFE","CNOT2","INO80E","NOTCH1","C17ORF49","BRD3","CACNG1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='D3PR1,D3PR2_vs_D3NR1,D3NR2 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2947.5951655440185,2806.238592778524,4339.2667369117335,3528.5613360203556),c(2627.5770785520353,3206.874312843782,4292.2522096658995,5205.314723141615),c(1124.6349914289704,2000.4957718437088,1567.1509081944462,2318.2610273594278),c(1955.5390958688702,1733.1071997465838,2671.1216590780896,3015.2474979456497),c(2211.553565462457,3139.8036007792857,3686.287191830714,4387.136632863044),c(3321.330574281085,2719.4938051751087,4530.8074034688325,4031.641215098819),c(1100.6336349045716,1006.954957128304,2894.0053437990773,3044.8721160926684),c(1637.806852355401,1483.6041508666576,1030.837041834569,977.6123988516174),c(2649.29259159792,3091.5126880928483,3421.612816224541,4162.528164366556),c(2333.8461915629646,2563.8897531854777,4440.2609065509305,3860.8956886878204),c(2101.8330784937766,3116.552420596927,3148.232046683954,3109.5076465952548),c(3091.6033046904113,2640.7975030194334,4341.008015698616,3843.659547220464),c(2102.976000233034,3198.8258273960423,2963.656495274386,3320.6503795703698),c(2259.5562785112543,1900.3368418273944,3978.8220280270107,3908.2950777230503),c(1650.3789914872289,2074.7206931950846,3815.141822060035,3380.976874706117),c(1410.3654262432412,731.5178995834393,1140.5376054081803,1475.8446131423866),c(1774.957461066251,2564.7840293463373,3322.359925372226,4585.352259737641),c(1076.6322783801727,1535.4721681965348,5067.121269828709,1967.6132743828973),c(2517.856591583355,1616.8512988347902,2545.749586422534,2720.0785753171726),c(4248.24010481858,2986.882377272234,7205.411620120687,4953.774783602383),c(3980.7964178324223,3181.8345803397033,6602.929159859266,6119.90747975321),c(1946.3957219548133,1703.5960864382052,3902.205761404171,1357.8847699751666),c(2302.9873046030234,1878.8742139667556,2772.115828717287,2490.6224420329913),c(2016.1139480494955,1740.2614090334632,2263.6624229475333,3607.2012314651693),c(994.3419131536627,2852.7409531432413,2655.450149996145,3570.035801426182),c(1770.3857741092224,1673.1906969689671,3141.2669315364233,2577.341778790628),c(960.0542609759502,990.8579862328249,1622.871829374693,2412.521176009033),c(678.8955131187076,1127.6822388443973,2098.240938193675,2114.120476855426),c(4967.137878811285,3286.4648911603176,8600.175928413744,5556.5011055390005),c(1756.6707132381375,1404.9078487109819,4835.531191173308,3291.564390844206),c(2851.5897394464237,3662.9551548823565,3691.511028191362,4309.035366839085),c(1518.942991472664,1333.365755842186,1769.1392474728414,2824.5726829630203),c(2496.1410785374705,3019.9705952240524,2770.374549930404,3384.7472806521014),c(2804.729948136883,3133.543667653266,3997.9760946827205,3428.376263741347),c(2046.9728350094367,2495.0304887992615,1582.8224172763908,3476.314282197432),c(2257.2704350327404,2713.2338720490893,2235.80196235741,3901.2928952519364),c(3617.3473047486696,4563.491248868327,5805.4234754669815,4781.41336892882),c(2927.0225742373914,1893.1826325405148,4024.095276485961,3215.6176425036674),c(690.3247305112785,1313.691680303267,2418.636234980095,2136.7429125313315),c(4291.671130910349,4805.840088461373,6320.841996384266,6717.786136902134),c(3946.50876565471,4217.406374615526,7240.237195858342,7259.647334282148),c(3379.6195829831963,2818.7584590305632,5707.911863401549,5467.627251097944),c(3960.2238265257947,4638.610446380562,5337.019481795531,6089.744232185337),c(2957.8814611973326,2676.568549453831,4917.371294156796,4088.7359337094367),c(1001.1994435892052,1615.0627465130704,2777.339665077935,2597.2710673622587),c(2883.591548145622,2716.810976692529,3303.205858716516,4319.269325835327),c(3435.6227482067934,4227.243412384985,6164.126905564822,5555.42384669729),c(2235.554921986856,1618.63985115651,3606.188367634109,3150.443482580226),c(717.7548522534485,1286.8633954774684,842.7789328512355,1518.9349668107775),c(1914.393913255615,2085.452007125404,4246.978961206949,4154.4487230537325),c(790.9018435659019,1343.2027936116453,1159.69167206389,1816.2584071226743),c(850.3337740072702,2191.870870267738,3952.70284622377,3401.983422119458),c(3369.3332873298828,4357.807731870538,3701.9587009126585,6620.832841148254),c(2258.4133567719973,2681.039930258131,3752.455785732257,3524.2523006535166),c(2356.7046263481066,3213.1342459698017,3957.926682584418,3524.2523006535166),c(1984.1121393502972,1371.8196307591638,2728.583859045219,2010.703628051288),c(2237.8407654653697,2772.256098665846,3789.022640256794,3184.3771360940837),c(1342.9330436270732,1587.340185526412,3200.4704102904357,4026.793550311125),c(970.340556629264,2093.5004925731437,1852.720629243212,1725.7686644190535),c(2410.4219480931893,3210.4514174872215,5444.978766582259,5278.568324377879))
targetgene="NON-TARGETING"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1548.6589566933483,1982.61024862651,6717.853559793526,7168.618962157673),c(3983.0822613109362,3501.0911697667057,8375.550964905873,10275.972091069509),c(2387.5635133080477,2622.0177036413743,7539.737147202169,9239.110455923854),c(950.9108870618936,1808.2263972588196,4330.56034297732,6785.114814508995),c(1786.3866784588217,1043.6202797235621,4509.91205802624,6732.86776068607),c(2390.992278525819,3114.763868275207,6932.0308505801,8079.979942244139))
targetgene="CBLL1"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1299.502017535304,1720.5873334945445,3726.3366039290163,3595.3513842063617),c(1701.8104697537974,1305.6431948555276,4565.632979206486,3613.6647845154275),c(3802.500626508317,3799.779407493929,9092.957825101554,7950.708881238967),c(3285.9000003641154,2947.534226194396,6178.057135859884,9629.616786043645),c(2857.304348142709,2771.361822504986,7466.603438153094,7502.569203087702),c(1323.5033740597028,1116.9509249140779,1899.7351564890453,3995.5530439015415),c(3266.470330796745,4103.833302186312,10759.361624148314,8346.601505567307),c(1944.1098784762992,1928.953678974913,5666.121172516364,6178.618086626393))
targetgene="SUZ12"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1722.383061060425,1667.8250400038073,2486.5461076685215,2308.56569778404),c(1736.09812193151,1876.1913854841757,2692.017004520682,3594.8127547855065),c(779.472626173331,1471.0842846146184,2068.639198816669,2772.325629140096),c(1320.0746088419316,1555.1462437354537,2288.0403259638915,2449.147976627165),c(1955.5390958688702,1791.2351502024806,2568.386210652009,3422.989969532798),c(795.4735305229302,1505.9610548881562,720.8894177694452,1448.9131420996423),c(2387.5635133080477,1775.1381793070013,1694.2642596368846,1900.8232261968915),c(2221.8398611157704,1553.3576914137336,1547.9968415387364,2081.264082183278),c(725.7553044282481,874.6020853210314,1607.2003202927488,1434.3701477365605),c(1296.073252317533,1570.3489384700727,3548.726167666979,3012.0157214205205),c(1461.79690450981,2019.275571221768,3355.4442223229976,3202.69053640315),c(1430.9380175498686,1958.4647922832912,3217.883198159263,4041.875174095062),c(873.1922087924119,810.214201739115,1702.9706535712983,2414.675693692452),c(1289.2157218819902,1056.1401459756014,1177.1044599327174,3272.712361114285),c(454.88285222431927,1104.4310586620386,1647.249732391051,2533.174166280527),c(2009.256417613953,1429.9475812150606,2979.3280043563304,3201.0746481405854),c(1635.5210088768868,1049.8802128495818,3341.5139920279357,1102.5744244899508),c(2088.1180176226917,2500.396145764421,4031.060391633492,2836.422530221828),c(465.169147877633,735.0950042268792,1253.720726555557,661.9755582306544),c(1365.791478412215,1598.9657756175911,651.2382662941366,2834.8066419592633),c(917.7661566234382,1083.8627069622598,1741.278786882718,1259.3155859587225),c(213.72636524107463,419.4155194433167,985.5637933756184,1486.6172015594843),c(2044.6869915309226,1243.043863595331,2718.1361863239226,1712.8415583185363),c(2537.2862611507257,2445.845299951964,3849.9673977976895,3640.057626137317),c(2401.2785741791326,2988.670929593954,5424.083421139667,4897.218694412621),c(1531.515130604492,1489.8640839926773,3491.2639676998497,1588.4181621010578),c(2006.970574135439,2213.333498128377,2810.4239620287067,2998.0113564782932),c(1596.661669742146,802.1657162913755,2991.5169558645093,1398.281976539283),c(2418.422400267989,2627.383360606534,4643.990524616209,3324.420785516354),c(621.7494261558535,703.7953385967809,741.7847632120379,758.390224563679),c(1445.7960001602107,1412.9563341587216,1481.828247637193,2232.080320022646),c(816.0461218295577,840.6195912083533,1500.9823142929029,713.1453532118686),c(723.469460949734,1396.8593632632424,1161.4329508507728,1399.359235380993),c(324.58977394901177,379.1730922046189,1168.3980659983038,401.2789185368898),c(385.16462612963716,1181.3388084959943,1547.9968415387364,1679.9851636463882),c(734.8986783423047,1123.2108580400975,851.4853267856491,948.5264101254536),c(938.3387479300657,959.5583206027267,2044.261295800311,1431.1383712114311),c(533.744452233058,1248.4095205604906,1250.2381689817914,1355.1916228708922),c(2880.162782927851,2564.7840293463373,3964.891797731949,3519.9432652866776),c(1278.9294262286767,570.5481906286483,959.4446115723777,491.2300318196557),c(1357.7910262374153,849.5623528169529,1074.369011506637,2964.0777029644355),c(920.0520001019523,1163.4532852787954,2155.703138160805,1203.8367556106693),c(1490.369947991237,1803.7550164545198,4240.013846059418,2797.1025824994213),c(1100.6336349045716,1197.4357793914735,1316.4067628833347,926.9812332912581),c(1250.3563827472494,969.3953583721861,1549.738120325619,1831.878660327466),c(1504.085008862322,1640.996755178009,4499.464385304943,1248.00436812077),c(1453.7964523350104,1497.018293279557,2817.3890771762376,2149.131389210994),c(1378.3636175440429,1650.8337929474683,823.6248661955256,1654.1309514453537),c(1061.7742957698306,943.4613497072476,1394.7643082930572,1786.6337889756555),c(316.5893217742122,1454.9873137191391,1053.4736660640444,1831.3400309066112),c(2544.143791586268,1306.5374710163874,2190.528713898459,2298.870368208652),c(1787.5296001980787,2071.143588551645,2711.1710711763917,2508.935842342057),c(2520.1424350618695,901.4303701468299,2204.458944193521,2692.608474853573),c(1916.6797567341293,1845.7859960149374,4051.955737076085,3652.9847322378346),c(1613.8054958310022,1344.9913459333652,2167.8920896689838,3214.0017542411024),c(896.0506435775536,1754.5698276072226,1405.2119810143533,1354.1143640291825),c(1670.9515827938562,1952.2048591572716,2052.9676897347244,2189.5285957751103),c(1641.235617573172,2831.2783252826025,4826.8247972388945,2944.68704381366),c(1734.955200192253,1167.9246660830952,1689.0404232762364,1665.4421692833064),c(885.7643479242398,1197.4357793914735,2230.5781259967616,1489.3103486637588))
targetgene="SAFE"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2254.984591554226,2197.236527232898,5185.528227336734,5084.123103449266),c(3666.4929395367244,1930.7422312966328,5178.563112189203,6170.000015892715),c(946.3392001048652,999.8007478414245,2080.828150324848,2872.510701419105),c(578.3184000640844,742.2492135137587,2822.612913536886,1642.2811041865464),c(613.748973981054,417.6269671215968,1373.8689628504644,1291.0947217891608),c(2330.4174263451937,2349.2634745790897,6843.225632449082,7509.032756137961))
targetgene="CNOT2"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1618.3771827880305,2202.602184198058,4095.4877067481525,5924.923629403742),c(2342.9895654770216,2514.7045643381803,5848.9554451390495,5879.678758051931),c(1642.3785393124292,1789.4465978807607,3602.7058100603435,3028.174604046167),c(906.3369392308673,2408.2857011958463,4307.923718747844,3827.5006645948174),c(1587.518295828089,1195.6472270697536,2106.947332128089,2648.440862343473),c(2297.272695906738,1918.2223650445935,5455.4264393035555,4749.634233098382),c(2737.297565520715,3375.8925072463126,6035.2722753355,6815.278062076868),c(1452.6535305957532,1165.2418376005153,2018.14211399707,2552.026196010448))
targetgene="INO80E"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1870.9628871638458,3668.3208118475163,4337.52545812485,3651.9074733961247),c(5164.863339702761,3823.0305876762877,7884.510347004947,6664.4618242375),c(3314.4730438455426,4230.820517028425,6468.850693269297,6352.595389562521),c(2737.297565520715,2289.346971801473,3935.2900583549426,4122.669587223295),c(2073.26003501235,2400.237215748107,3237.037264814973,3161.2160709973236),c(4594.545391813476,4079.6878458430933,7365.609268513897,5631.909224458685),c(860.620069660584,1453.1987613974193,4689.263773075159,2666.7542626525387),c(1394.364521893642,1649.9395167866082,3479.0750161916703,3235.008301654443),c(2248.1270611186833,2034.478265956387,2115.6537260625023,2738.3919756262385),c(3567.058748221358,4564.385525029186,5100.205566779481,6292.268894426774),c(4444.822643970798,4311.30537150582,7431.777862415441,8650.388498929464),c(1542.944347997063,1484.4984270275174,1854.4619080300947,2380.7420401785944))
targetgene="NOTCH1"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3988.796870007222,2808.921421261104,4750.208530616054,5586.125723686019),c(3146.4635481747514,5343.300061138203,5936.019384483186,6471.093862150596),c(4769.41241791981,6146.360053590439,7835.754540972231,10195.177677941276),c(2591.0035828958084,3062.00157478447,4583.045767075314,5036.185084993181),c(2584.146052460266,2892.9833803819392,4518.618451960653,4097.354004443115),c(1986.3979828288113,2819.6527351914233,3883.051694748461,4344.046279194653))
targetgene="C17ORF49"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1906.3934610808155,1382.5509446894832,3623.6011555029363,2363.505898711238),c(1865.2482784675606,2354.629131544249,2972.3628892087995,3627.669149457655),c(5135.147374482077,4724.4609578231175,8392.9637527747,6626.219135356803),c(3663.064174318953,3814.087826067688,6028.307160187969,6887.454404471423),c(2946.4522438047616,2786.5645172396053,4708.41783973087,4355.89612645346),c(1987.5409045680685,2413.651358161006,5218.612524287506,3132.6687116920148))
targetgene="BRD3"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3085.888695994126,2706.07966276221,4978.31605169769,5243.557412022312),c(1507.5137740800933,1122.3165818792377,1803.9648232104958,2557.9511196398516),c(4613.9750613808465,5070.545832075918,7069.591874743835,8007.803599849584),c(2828.731304661282,3538.6507685228235,5133.289863730252,6029.41773704959),c(2482.4260176663856,3278.4164057125777,4950.455591107567,4042.4138035159167),c(5257.440000582585,5260.132378178228,11025.77727854137,10782.822376093956))
targetgene="CACNG1"
collabel=c("D3NR1","D3NR2","D3PR1","D3PR2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("sample1_summary.Rnw");
library(tools);

texi2dvi("sample1_summary.tex",pdf=TRUE);

