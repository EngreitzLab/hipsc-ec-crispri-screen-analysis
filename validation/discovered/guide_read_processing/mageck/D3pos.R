pdf(file='D3pos.pdf',width=4.5,height=4.5);
gstable=read.table('D3pos.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("PLCG1","SMAD1","EIF3M","SETDB1","HSD17B12","ASCC3","PSMD4","UPF2","NAA25","SRPR")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='D3PR1,D3PR2_vs_D0R1,D0R2 neg.'


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
targetmat=list(c(2233.0771257545407,2554.4494819603347,477.7270803828376,93.12293389547156),c(1518.5613145069426,811.6688577110946,326.2798008209124,63.492909474185154),c(2708.2731833553526,2240.9118462893266,698.2165020979934,302.2967967743149),c(2363.9282140794016,3281.9653811796434,140.31145018237189,145.68095340465817),c(1942.1056267163622,1957.234937825081,288.4179809304311,273.72498751093156),c(2456.9013557839085,1516.1105586342687,347.4378766420637,601.7716864608882))
targetgene="PLCG1"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(1530.613388431601,2334.565945255991,135.85711843055054,816.9421019011824),c(1356.7191789472458,1517.4678644163944,511.13456852149756,427.5189237928467),c(2362.206489233022,2429.5773500047817,412.025687043473,733.6958428128062),c(905.6272691957504,1411.5980134105994,91.31380091233726,411.6456964243004),c(1327.44985655879,1924.6595990540673,550.1099713499342,431.75178442445906),c(2744.4294051293273,1150.9953032424887,384.18611359458964,667.3810262508795))
targetgene="SMAD1"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(764.4458317926106,1692.5603103105936,202.67209470787049,292.7728603531871),c(2190.0340045950466,1742.78062424924,658.1275163316014,396.12520744172184),c(1160.442546459954,2000.6687228530996,90.20021797438193,214.46493866835874),c(655.9771664706861,844.2441964821085,18.93090994524065,165.4343030188491),c(786.8282547955474,1032.9097001975636,134.7435354925952,284.6598774759301),c(752.3937578679523,1688.488392964217,290.64514680634176,503.71041516186887))
targetgene="EIF3M"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(578.4995483835972,803.5250230183411,75.7236397809626,188.0095597207816),c(1597.7606574404115,1981.6664419033414,139.19786724441656,220.81422961577726),c(1360.1626286400053,842.8868906999829,247.21541222608377,138.97892407127196),c(869.4710474217756,1421.0991538854785,131.40278667872923,441.62845923155453),c(1196.5987682339287,1516.1105586342687,222.71658759106649,170.37264042239684),c(1678.6817252202598,1161.8537494994932,563.4729666053981,44.7977750178973))
targetgene="SETDB1"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(1735.4986451507916,2274.844490842466,468.8184168791949,338.62885052898747),c(1525.4482138924618,3349.830670285922,616.9249476272541,553.0937891973463),c(451.09190975149545,709.8709240516764,23.38524169706198,163.31787270304292),c(1735.4986451507916,1750.9244589419934,484.40857801056956,883.6096568490767),c(1628.7517046752469,1376.3080630753343,589.0853741783708,378.84102652930477),c(955.5572897407633,1250.0786253376557,253.89690985381577,211.99576996658487))
targetgene="HSD17B12"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(3255.7816845041143,2345.424391512996,1054.5630422436998,1177.7934707461345),c(1200.0422179266882,1398.0249555893436,153.67444543783586,372.13899719591853),c(1525.4482138924618,1802.5020786627654,881.9576868606232,621.1722976891115),c(2208.972977905224,2080.7497639985086,545.6556395981129,719.9390457600662),c(1336.058480790689,2287.060242881596,515.5889002733189,341.09801923076134),c(2444.8492818592504,2747.1869030221665,155.90161131374654,719.5863073740984))
targetgene="ASCC3"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(919.4010679667884,1317.9439144439345,40.088985766391964,129.1022492641765),c(860.8624231898768,935.1836838845221,148.1065307480592,185.18765263304005),c(655.9771664706861,602.6437672637559,31.180322262749307,225.0470902473896),c(328.84944565853294,396.33328838066825,37.8618198904813,8.112982877256991),c(612.9340453111923,475.0570237439517,40.088985766391964,45.855990175800386),c(576.7778235372175,625.7179655598907,6.681497627731995,8.818459649192382))
targetgene="PSMD4"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(2472.3968794013263,3761.094322269972,883.0712697985786,1196.4886052024224),c(1060.5825053699282,1464.5329389134968,37.8618198904813,190.47872842255546),c(761.0023820998512,758.7339322081972,59.019895711632614,8.112982877256991),c(1349.8322795617269,1864.9381446405419,181.51401888671919,279.72154007238237),c(1633.9168792143862,1938.2326568753228,396.43552591209834,384.13210231882016),c(1036.4783575206118,1164.5683610637445,1015.5876394152631,934.7567228143926))
targetgene="UPF2"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(2291.615770531452,2114.6824085516478,761.6907295614474,518.5254273725121),c(5686.857167592327,4001.337445706199,1530.0629567506267,1399.3131771338471),c(1382.545051642942,970.4736342197871,523.3839808390062,406.0018822488173),c(2537.8224235637567,3481.489331152103,613.5841988133882,685.370683935232),c(680.0813143200027,2373.927812937633,677.058426276842,268.7866501073838),c(2594.639343494289,2344.06708573087,870.8218574810699,940.7532753758434))
targetgene="NAA25"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(690.4116633982812,1750.9244589419934,101.33604735393524,195.06432744013551),c(2310.5547438416293,1855.4370041656628,773.940141878956,328.752175721892),c(1048.5304314452699,3018.648059447282,624.7200281929414,240.5675792299682),c(835.0365504941806,540.2077012859794,1203.7831559297142,585.5457207063743),c(1406.6491994922587,560.567288017863,106.90396204371191,116.05092898337175),c(1747.55071907545,1370.8788399468322,408.684938229607,449.7414421088115))
targetgene="SRPR"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetgenelist=c("NON-TARGETING","CBLL1","SAFE","APAF1","ACVR1B","SOX4","NOTCH1","SUZ12","C17ORF49","BAX")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='D3PR1,D3PR2_vs_D0R1,D0R2 pos.'


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
targetmat=list(c(1194.877043387549,1121.134576035726,2775.0486813846883,2310.789166474372),c(855.6972486507376,1506.6094181593896,2744.981942059894,3408.863761991807),c(1132.8949489178779,913.4667913705129,1002.2246441597991,1518.1860132049605),c(2040.243942960008,1275.8674351980417,1708.23622682348,1974.6294846471583),c(1618.4213555969684,1178.1414188850001,2357.4550796514386,2873.0541537068784),c(2296.7809450705913,974.5455515661638,2897.5428045597746,2640.2468189681995),c(494.1350309109893,650.149469638151,1850.7748428817624,1994.0300958753814),c(814.3758523376234,1303.0135508405533,659.2410992695568,640.220170531367),c(1206.9291173122074,1556.829732098036,2188.190473082228,2725.9622467583495),c(1645.9689531390445,1923.3022932719416,2839.6364917860974,2528.42875061644),c(1613.2561810578293,818.4553866217225,2013.357951823241,2036.358702191505),c(1126.0080495323589,1074.9861794434564,2776.1622643226438,2517.1411222654738),c(607.768870772053,1687.1310871820913,1895.3181603999758,2174.6321494908416),c(1074.3563041409664,1167.2829726279956,2544.5370132279345,2559.469728581597),c(3030.2357296283662,2125.5408548086525,2439.860217060133,2214.1388487192235),c(785.1065299491677,1227.004427041521,729.3968243607427,966.5031775514851),c(989.9917866683584,1560.9016494444127,2124.7162456187743,3002.86187974299),c(1112.234250761321,1505.2521123772642,3240.5263494500173,1288.5533239399908),c(919.4010679667884,822.5273039680992,1628.058255290696,1781.3288491368612),c(1492.7354418112463,2368.4985898091304,4608.006197259166,3244.1349357448935),c(1561.6044356664365,2177.1184745294245,4222.70650072662,4007.813541364954),c(1096.7387271439031,707.1563124874252,2495.5393639578997,889.2534710245599),c(900.4620946566112,879.5341468173734,1772.824037224889,1631.062296714623),c(638.7599180068886,677.2955852806625,1447.657819341932,2362.2889708256553),c(671.472690088104,849.6734196106107,1698.2139803818818,2337.9500221938843),c(1048.5304314452699,1564.9735667907894,2008.9036200714195,1687.853176855422),c(284.08459965265934,799.4531056719644,1037.8592981743698,1579.9152307493073),c(664.5857907025849,804.8823288004668,1341.8674402361755,1384.498164923204),c(1923.1666534061849,2436.3638789154093,5499.986130561387,3638.849189642745),c(928.0096921986872,667.7944448057835,3092.419818701958,2155.584276648586),c(1511.6744151214236,1843.2212521265326,2360.7958284653046,2821.9070877415625),c(1470.3530188083096,598.5718499173792,1131.4002649626177,1849.7600960145942),c(1101.9039016830422,1160.4964437173678,1771.7104542869338,2216.6080174209974),c(2200.3643536733252,2219.194953775317,2556.7864255454433,2245.1798266843807),c(1282.6850105529165,602.6437672637559,1012.2468906013971,2276.5735430355053),c(519.9609036066856,1566.330872572915,1429.8404923346468,2554.8841295640173),c(1962.766324872919,2633.173217323618,3712.685515143078,3131.2586522352312),c(933.1748667378265,2075.320540870006,2573.490169614773,2105.848164227141),c(1129.4514992251184,447.91090810144016,1546.7667008199567,1399.3131771338471),c(3198.964764573582,1879.868508243923,4042.3060647778566,4399.353149789095),c(1849.1324850118554,1711.5625912603516,4630.277856018272,4754.207966072598),c(1733.776920304412,2155.401582015415,3650.3248706175796,3580.647355958075),c(2202.086078519705,1533.7555338019013,3413.131704833094,3988.060191750763),c(2165.92985674573,1338.3035011758182,3144.7582167858586,2677.637087880775),c(692.133388244661,646.0775522917743,1776.1647860387552,1700.9044971362268),c(1546.1089120490187,1183.5706420135025,2112.4668333012655,2828.6091170749487),c(926.2879673523074,1704.7760623497238,3942.0836003618765,3638.1437128708094),c(1187.99014400203,1438.7441290531108,2306.2302645054933,2063.16681952505),c(533.7347023777236,449.26821388356575,538.9741419703809,994.7222484289008),c(1587.4303083621328,897.179121985006,2716.0287856730556,2720.671170968834),c(297.8583984236974,594.4999325710024,741.6462366782514,1189.4338374830686),c(588.8298974618757,701.727089358923,2527.8332691586043,2227.8956457719637),c(2019.583244803451,2337.2805568202425,2367.4773260930365,4335.860240314911),c(1353.2757292544864,861.8891716497409,2399.771231293741,2307.9672593866303),c(1473.796468501069,989.4759151695451,2531.1740179724707,2307.9672593866303),c(426.98776190217893,1112.9907413429726,1744.984463776006,1316.7723948174066),c(798.8803287202057,665.0798332415324,2423.156472990803,2085.389337841015),c(810.9324026448639,1020.6939481584334,2046.765439961901,2637.07217349449),c(984.8266121292191,779.0935189400808,1184.8522459844737,1130.1737886404958),c(1150.1121973816755,1145.5660801139863,3482.1738469863244,3456.8361824834137))
targetgene="NON-TARGETING"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(2326.050267459047,895.8218162028803,4296.202974631672,4694.595178844057),c(946.9486655088645,1748.2098473777423,5356.333931565148,6729.5429274916905),c(1825.028337162539,1612.4792691651846,4821.814121346589,6050.521534503878),c(468.309158215293,682.7248084091649,2769.4807666949114,4443.445448035058),c(607.768870772053,878.1768410352479,2884.1798093043108,4409.229824596191),c(1348.110554715347,1258.2224600304091,4433.173676000178,5291.4285279013975))
targetgene="CBLL1"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(552.6736756879009,1715.6345086067283,1590.1964354002146,1511.836722257542),c(1126.0080495323589,1017.9793365941822,1721.5992220789437,2354.1759879483984),c(432.15293644131816,718.0147587444299,1322.9365302909348,1815.5444725757277),c(1356.7191789472458,722.0866760908066,1463.2479804733068,1603.9014409951105),c(676.6378646272432,789.9519651970854,1642.5348334841153,2241.6524428247035),c(626.7078440822303,678.6528910627882,461.02333631350757,948.8662582531003),c(1573.6565095910948,1069.556956314954,1083.5161986305384,1244.8137640799966),c(735.1765094041548,1541.8993684946547,989.9752318422904,1362.9811233791747),c(309.91047234835565,78.72373536328342,1027.8370517327717,939.3423218319726),c(755.8372075607118,792.6665767613365,2269.4820275529673,1972.5130543313521),c(1019.2611090568141,498.13122204008647,2145.8743214399256,2097.382442963916),c(977.9397127437,720.729370308681,2057.9012693414543,2646.9488483015857),c(974.4962630509406,559.2099822357375,1089.084113320315,1581.3261842931781),c(304.74529780921637,363.7579496096544,752.7820660578046,2143.2384331397166),c(168.72903494521586,1146.923385896112,1053.4494593057443,1658.928629206071),c(1129.4514992251184,880.891452599499,1905.3404068415737,2096.3242278060134),c(203.16353187281092,986.761303605294,2136.9656579362827,722.0554760758723),c(1279.241560860157,944.6848243594011,2577.9445013665945,1857.5203405058835),c(89.52969201174719,260.60271016811066,801.7797153278393,433.5154763542975),c(535.4564272241034,604.0010730458814,416.4800187952943,1856.4621253479804),c(387.38809043544455,491.3446931294586,1113.5829379553325,824.7023463924717),c(588.8298974618757,177.8070574584505,630.2879428827181,973.557945270839),c(1148.3904725352957,586.356097878249,1738.3029661482738,1121.7080673772712),c(2105.6694871224386,883.6060641637501,2462.1318758192397,2383.806012369685),c(1179.3815197701313,1259.5797658125348,3468.81085173086,3207.0974052182855),c(890.1317455783326,423.4794040231798,2232.7337906004414,1040.2255002187335),c(633.5947434677494,757.3766264260715,1797.3228618599064,1963.341856296192),c(395.9967146673433,249.74426391110603,1913.135487407261,915.708849972137),c(1181.103244616511,1521.539781762771,2969.9256955268716,2177.1013181926155),c(344.34496927595075,435.69515606231,474.38633156897157,496.655647442515),c(561.2822999197997,563.2818995821142,947.6590801999879,1461.7478714501294),c(533.7347023777236,299.96457784975235,959.9084925174965,467.02562302122857),c(463.1439836761537,411.2636519840496,742.7598196162066,916.4143267440724),c(91.25141685812694,92.29679318453918,747.2141513680281,262.790097545933),c(153.23351132779808,321.6814703637616,989.9752318422904,1100.1910258332416),c(483.8046818327108,424.83670980530536,544.5420566601575,621.1722976891115),c(590.5516223082554,897.179121985006,1307.3463691595603,937.2258915161664),c(874.6362219609149,458.7693543584448,799.5525494519286,887.4897790947214),c(1222.4246409296252,1441.458740617362,2535.6283497242916,2305.1453522988886),c(866.0275977290161,367.8298669560312,613.5841988133882,321.6974080025381),c(557.8388502270402,962.3297995270336,687.08067271844,1941.1193379802273),c(724.8461603258763,251.10156969323162,1378.6156771887015,788.370292637799),c(1153.555647074435,663.7225274594068,2711.574453921234,1831.7704383302416),c(609.4905956184327,1119.7772702536004,841.8687010942313,607.0627622504036),c(817.819302030383,544.2796186323561,991.0888147802458,1199.6632506761316),c(817.819302030383,844.2441964821085,2877.498311676579,817.29484028715),c(680.0813143200027,308.1084125425058,1801.7771936117279,1407.4261600111042),c(1020.9828339031939,413.9782635483008,526.7247296528722,1083.2595833067924),c(246.20665303230479,814.3834692753458,891.9799333022212,1170.0332262548452),c(275.47597542076056,153.37555338019013,673.7176774629761,1199.310512290164),c(1095.0170022975233,965.0444110912847,1400.8873359478082,1505.4874313101236),c(1079.5214786801055,811.6688577110946,1733.8486343964526,1643.0554018375246),c(499.30020545012854,1258.2224600304091,1409.7959994514508,1763.3391914525089),c(1239.6418893934226,1653.1984426289519,2591.3074966220584,2392.2717336329097),c(1048.5304314452699,582.2841805318723,1386.410757754389,2104.789949069238),c(401.1618892064826,584.9987920961233,898.6614309299532,886.784302322786),c(1358.4409037936257,1026.1231712869358,1312.914283849337,1433.8815389586814),c(1680.4034500666396,1104.8469066502191,3086.8519040121814,1928.4207560853902),c(1060.5825053699282,650.149469638151,1080.1754498166724,1090.667089412114),c(1067.4694047554472,422.12209824105423,1426.4997435207808,975.3216372006775))
targetgene="SAFE"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(1694.1772488376776,1275.8674351980417,2940.972539140033,3200.042637498932),c(1127.7297743787387,1210.7167576560141,3791.7499037379066,2299.5015381234057),c(1162.1642713063336,1577.1893188299196,3237.1856006361513,3353.1310970089116),c(1713.1162221478548,765.520461118825,2762.7992690671795,3777.8281137140166),c(1048.5304314452699,985.4039978231684,3189.301534304072,2030.362149630054),c(859.1406983434971,1744.1379300313656,2589.080330746148,3069.5294346908845))
targetgene="APAF1"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(2245.129199679199,3381.0487032748106,6606.887570888987,5912.600825590509),c(1143.2252979961563,1450.959881092241,3062.353079377164,2866.352124373492),c(1463.4661194227906,1270.4382120695393,2244.98320291795,3995.4676978560847),c(793.7151541810664,1007.1208903371776,2673.712634030753,2127.717944157138),c(1382.545051642942,1498.4655834666362,2758.344937315358,3029.669997076535),c(370.170841971647,606.7156846101326,1761.688207845336,1982.7424675244154))
targetgene="ACVR1B"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(1931.7752776380835,2008.812557545853,3193.7558660558934,4074.1283579268807),c(1596.0389325940316,1400.7395671535946,3865.2463776429586,3250.4842266923124),c(1053.6956059844092,1112.9907413429726,2108.012501549444,2437.7749854227422),c(1119.12115014684,489.987387347333,2023.380198264839,1747.4659640839625),c(2038.5222181136282,2019.6710038028575,4699.319998171502,4439.565325789413),c(1439.361971573474,1041.053534890317,3023.3776765487273,3242.723982201023))
targetgene="SOX4"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(988.2700618219786,722.0866760908066,2773.935098446733,2391.566256860974),c(2496.501027250643,2473.0111350328,5042.303543061745,4364.432049578294),c(2165.92985674573,1294.8697161477999,4136.96061450406,4160.196524102998),c(1787.1503905421844,1188.9998651420049,2516.6974397790514,2699.8596061967396),c(957.279014587143,955.5432706164057,2070.1506816589626,2070.2215872444035),c(2410.414784931655,2630.458605759367,4710.455827551056,3688.232563678222),c(513.0740042211665,899.893733549257,2998.87885191371,1746.4077489260594),c(1324.0064068660306,893.1072046386292,2224.938710034754,2118.546746121978),c(986.5483369755989,2535.4472010105765,1353.0032696157289,1793.3219542597628),c(2651.4562634248205,2273.4871850603404,3261.6844252711685,4120.6898248746165),c(2177.9819306703885,1828.2908885231514,4752.771979193359,5664.9784786411865),c(812.6541274912437,971.8309400019126,1185.965828922429,1559.1036659772133))
targetgene="NOTCH1"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(879.8013965000541,1662.699583103831,2383.0674872244113,2354.528726334366),c(1348.110554715347,1501.1801950308873,2919.8144633188813,2366.5218314572676),c(2419.023409163554,2858.485977156464,5815.130102002746,5206.77131526915),c(3181.7475161097846,3412.266736263699,3950.9922638655194,6306.256864330457),c(1694.1772488376776,2035.9586731883644,4775.043637952465,4913.292978144028),c(1776.8200414639057,864.6037832139921,1214.9189853092676,2616.6133471083635),c(3131.817495564772,2672.5350850052596,6880.8289736259985,5466.034028955406),c(1194.877043387549,1491.6790545560084,3623.5988801066514,4046.262025435433))
targetgene="SUZ12"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(2338.1023413837056,1769.9267398917514,3037.8542547421466,3658.249800870968),c(1978.261848490337,2384.7862591946378,3796.204235489728,4237.798969015891),c(2432.797207934592,2857.128671374338,5011.123220798996,6676.632169596537),c(1484.1268175793477,1645.0546079361984,2930.9502926984346,3298.103908797951),c(793.7151541810664,2219.194953775317,2889.7477239940877,2683.280902056258),c(1957.60115033378,1114.3480471250982,2483.2899516403913,2844.8350828294624))
targetgene="C17ORF49"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
targetmat=list(c(979.6614375900798,1159.1391379352422,2096.8766721698908,2869.526769847201),c(3217.9037378837597,4245.652486488802,4934.285998080078,4828.2830271258135),c(1243.085339086182,941.9702127951499,3315.1364062930243,2469.1687017738673),c(2785.7508014424416,4009.4812803989525,7299.536158297204,8901.353169894792),c(533.7347023777236,785.8800478507087,1846.320511129941,2796.509923951888),c(2506.831376328921,2297.9186891386007,3642.529790051892,5233.226694216727))
targetgene="BAX"
collabel=c("D0R1","D0R2","D3PR1","D3PR2")

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
Sweave("D3pos_summary.Rnw");
library(tools);

texi2dvi("D3pos_summary.tex",pdf=TRUE);

