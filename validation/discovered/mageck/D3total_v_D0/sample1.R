pdf(file='sample1.pdf',width=4.5,height=4.5);
gstable=read.table('sample1.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("SSRP1","ZC3H8","GOLT1B","UPF2","EIF3M","SMC2","COPS6","EIF4E","HSD17B12","SETDB1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='D3TOTAL1,D3TOTAL2_vs_D0R1,D0R2 neg.'


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
targetmat=list(c(2472.537316315824,3586.0437607089125,690.7378219278561,936.5343360279592),c(3325.0475785038834,4377.530499065718,523.6335628239119,984.7513596453178),c(947.8059712543991,1714.8879331064124,345.06727872727834,676.8131904080767),c(3232.3274291420403,4010.199474341149,337.22484057438567,524.1752506132482),c(4486.625005231422,2546.9637349686955,297.4093853366228,420.9375742791491),c(2861.446831694667,1757.5064497871635,2523.455291811548,1938.1468634413104))
targetgene="SSRP1"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1671.5382482176765,1286.6733131236278,194.25116040241892,287.5272819390955),c(2279.37033847865,1966.5401268403712,545.3510838626917,435.4322623604409),c(1594.271457082807,1037.050572564943,241.9090537930745,152.93374975567116),c(937.5037324364164,1215.6424519890427,127.89206833948079,106.78739586421139),c(1310.9598895882857,1158.8177630813746,196.66421829561668,452.589240089317),c(4429.962691732517,3348.5977392018704,213.55562354800094,374.1996004660039))
targetgene="ZC3H8"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(2191.8013085257976,2695.1138167636873,96.52231572791003,201.15077337302978),c(602.680970851982,1377.9987060109515,506.1388930982282,494.5942545289791),c(1486.0979494939897,1680.387229126757,270.8657485114475,286.6398520565674),c(3979.2397434457785,5745.381939200301,749.857740311201,640.4285652244257),c(1918.791979849259,1611.3858211674456,168.3107880505431,226.5904300055012),c(1099.7639938196423,2725.555614392795,456.0679418143749,417.97947467072214))
targetgene="GOLT1B"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(3698.5037356557527,5623.614748683869,1360.361387290232,1563.6514530144639),c(1586.5447779693202,2189.779976120496,492.2638102123411,271.2577340927475),c(1138.397389387077,1134.4643249780884,139.95735780546954,32.243285731853305),c(2019.2388083245892,2788.468662826285,226.82744196058857,331.60296610465645),c(2444.2061595663713,2898.0591342910734,612.3134403989292,457.9138193844854),c(1550.486942106381,1741.2708243849727,898.2608007428627,1120.2323217112703))
targetgene="UPF2"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1143.5485087960683,2530.7281095665044,444.6059168216856,415.02137506229525),c(3276.111944118466,2605.8178770516374,738.9989797918112,539.8531785379108),c(1735.9272408300678,2991.413980353671,271.46901298474694,316.2208481408365),c(981.2882474128425,1262.3198750203414,85.0602907352207,211.20831204168127),c(1177.0307849545118,1544.4138663834083,492.2638102123411,382.48227936959927),c(1125.5195908645987,2524.639750040683,517.6009180909175,570.0257945438652))
targetgene="EIF3M"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(2992.800376623945,868.6059590172125,59.11991838334489,341.06888485162256),c(1921.3675395537546,2187.750522945222,611.7101759256298,226.88623996634388),c(865.3880607105382,994.4320558841919,111.60392756039597,115.36588472864943),c(839.6324636655818,1136.4937781533622,183.39239988302904,277.1739333096013),c(1756.531718466033,1723.005745807508,225.01764854069026,201.74239329471516),c(1697.293845262633,1219.7013583395903,142.97368017196672,194.34714427364787))
targetgene="SMC2"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1980.6054127571545,2230.369039625973,482.61157863955015,382.77808933044196),c(811.3013069161296,2047.7182538513257,437.36674314209233,183.9937956441537),c(2135.1389950268936,1398.2932377636903,333.60525373458904,337.8149752823529),c(2730.0932867653887,2595.6706111752683,895.2444783763655,723.2553542603791),c(2809.9356376047535,3060.4153883129825,647.906044323596,651.0777238147625),c(1478.371270380503,1442.9412076197152,364.37174187286035,649.0070540888637))
targetgene="COPS6"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1450.0401136310506,2352.1362301424047,941.0925783471228,794.2497448626249),c(2196.952427934789,2711.3494421658784,786.6568731824667,619.7218679654374),c(2377.2416072494843,2546.9637349686955,818.0266257940375,306.7549293938704),c(1426.8600762905899,2287.193728533641,648.5093087968955,433.0657826736994),c(1988.3320918706415,2181.6621634194007,243.7188472129728,587.7743921944267),c(2789.3311599687886,671.749001015648,43.435042077559515,69.51534079803235))
targetgene="EIF4E"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(2596.164182131615,3401.363521758991,652.7321601099916,589.5492519594828),c(2281.9458981831453,5008.690436575889,834.9180310464217,701.9570370797054),c(674.7966425778602,1061.4040106682294,202.69686302861106,257.94628585482644),c(2596.164182131615,2617.9945961032804,1388.1115530620061,1207.7920701207067),c(2436.4794804528847,2057.865519727695,970.6525375387952,722.9595442995364),c(1429.4356359950855,1869.126374427226,714.8684008598336,449.3353305200474))
targetgene="HSD17B12"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(865.3880607105382,1201.4362797621257,238.8927314265773,442.2358914598228),c(2390.1194057719626,2963.001635899837,621.9656719717202,565.8844550920676),c(2034.692166551563,1260.2904218450676,901.2771231093599,347.28089402931903),c(1300.657650770303,2124.8374745117326,214.15888802130038,767.6268483867827),c(1790.0139946244765,2266.8991967809025,710.0422850734382,726.5092638296487),c(2511.1707118832583,1737.211918034425,1033.3920427619366,237.83120851752344))
targetgene="SETDB1"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetgenelist=c("APAF1","CASP3","NON-TARGETING","BAX","DBR1","TBX3","CBLL1","SMAD5","CYCS","SAFE")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='D3TOTAL1,D3TOTAL2_vs_D0R1,D0R2 pos.'


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
targetmat=list(c(2534.350749223719,1907.6859847574292,5378.102779464487,5478.4004748066345),c(1686.9916064446504,1810.272232344284,5671.289313488013,4446.023711465644),c(1738.5028005345634,2358.2245896682266,6623.8439168278255,5286.715620180571),c(2562.681905973171,1144.6115908544575,4634.27768388628,5565.9602232160705),c(1568.5158600378506,1473.383005248823,5413.092118915854,4093.1224281803134),c(1285.204292543329,2607.8473302269113,4652.978882558563,4873.764914844174))
targetgene="APAF1"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(718.5811575542862,1763.5948093129853,3254.6118334504663,4014.7327885570003),c(2124.836756208911,3096.945545467912,5741.871256864048,5126.978241325518),c(2271.643659365163,2642.348034206567,7579.414842534135,8135.365543095683),c(2681.1576523799713,3411.5107876353604,5927.676714640274,7318.338431248171),c(1596.8470167873027,2301.3999007605585,5390.771333403775,4958.070753684341),c(607.8320902609732,819.8990828106398,3123.7034427444883,2536.8662241869165))
targetgene="CASP3"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1787.4384349199809,1676.3283227762092,3059.1541441014483,2866.1027106048314),c(1280.0531731343378,2252.6930245539857,2873.951950798521,3919.481981165654),c(1694.7182855581375,1365.8219869593083,1136.5502676961405,1934.8929538720408),c(3052.0382498273448,1907.6859847574292,1957.5932158566752,2229.2238649105184),c(2421.0261222259105,1761.5653561377112,2444.4276458093213,3447.9609035824046),c(3435.7966457971966,1457.1473798466322,3322.7807189333025,3113.695647830164),c(739.1856351902514,972.1080709561794,1583.5692424110239,2005.2957245526013),c(1218.2397402264423,1948.2750482629065,1221.6105584313614,1027.643803967508),c(1805.4673528514504,2327.7827920391187,2583.7817391414915,3308.634412025497),c(2462.235077497841,2875.735149363061,2770.190461391018,2968.4529570564027),c(2413.2994431124234,1223.7602646901382,2200.105534123049,2738.608617481632),c(1684.4160467401548,1607.326914816898,3135.768732210477,2984.426694941908),c(909.1725756869643,2522.610296865409,2136.762764426608,2881.7806385294944),c(1607.1492556052854,1745.3297307355203,2571.113185202203,2774.9932426652827),c(4532.985079912343,3178.1236724788664,2192.866360443456,2543.078233364613),c(1174.4552252500162,1834.6256704475702,1139.5665900626377,1052.491840678294),c(1480.9468300849985,2333.87115156494,2087.8983420893537,3366.6131643506646),c(1663.8115691041896,2250.663571378712,2323.774751149434,1588.4994897252498),c(1375.3488822006768,1229.8486242159597,2210.964294642439,2028.6647114591738),c(2233.010263797728,3541.3957908528873,4738.642437767083,3708.569479084815),c(2336.032651977554,3255.242893139273,4388.749043253409,4413.484615772947),c(1640.6315317637288,1057.3451043176815,2379.2750826929823,1309.2548866897496),c(1347.0177254512248,1315.085657577462,2175.974955191072,1989.3219866670959),c(955.5326503678859,1012.6971344616566,1848.402346189477,2556.685491563377),c(1004.4682847533034,1270.437687721437,1444.8184135521533,2904.2621955535387),c(1568.5158600378506,2339.9595110907617,2022.7457789730145,1968.9110993689503),c(424.9673512417822,1195.3479202363042,1068.9846466866036,1652.6902512281138),c(994.1660459353208,1203.4657329373995,1085.2727874656885,1534.0704569301947),c(2876.9001899216405,3642.86844961658,5601.310634585278,4138.6771621500875),c(1388.2266807231551,998.4909622347395,2602.4829378137742,2272.412119193551),c(2261.3414205471804,2755.997412021903,2784.065544276905,3578.1172863531883),c(2199.5279876392847,894.9888502957726,1414.6551898871812,1992.2800862755228),c(1648.3582108772157,1735.1824648591512,2277.3233867053773,2857.820031701236),c(3291.56530234544,3318.1559415727625,2865.506248172329,2919.348503556516),c(1918.791979849259,901.0772098215942,1628.8140779084817,2734.4672780298342),c(777.8190307576862,2341.9889642660355,1966.0389184828673,3040.038967580334),c(2936.1380631250404,3937.1391600312904,3920.6158119730453,4135.423252580818),c(1395.953359836642,3103.0339049937334,2939.1045139148605,2392.215153334841),c(1689.567166149146,669.7195478403742,1202.3060952857793,1608.0229471408675),c(4785.389930952917,2810.7926477542974,4455.108135316347,5279.024561198661),c(2766.151122628328,2559.1404540203384,4591.44590628202,5381.966427571917),c(2593.588622427119,3222.771642334891,3761.353991021994,3935.1599090903164),c(3294.1408620499355,2293.2820880594627,3939.317010645328,4878.7936841785),c(3240.054108255527,2001.040830820027,3264.8673294965565,3130.85262555904),c(1035.3750012072512,966.0197114303578,1490.6665135229105,1960.628420465355),c(2312.852614637093,1769.6831688388068,2666.4289719835147,3270.7707370376324),c(1385.6511210186595,2548.9931881439693,3948.969242218119,4449.277621034913),c(1777.1361961019982,2151.2203657902924,2429.3460339768353,2265.608490094169),c(798.4235083936513,671.749001015648,670.8300943089747,1259.8546232290203),c(2374.666047544989,1341.4685488560222,2481.8300431538864,2971.4110566648296),c(445.57182887774735,888.900490769951,819.2331547406363,1441.7777491472752),c(880.8414189375121,1049.2272916165862,1818.239122524505,2593.3659267078706),c(3021.131533373397,3494.7183678215883,3060.963937521347,5077.577977864788),c(2024.3899277335806,1288.7027662989017,2492.085539199977,2822.3228364001134),c(2204.679107048276,1479.4713647746446,2615.151491753062,2998.3297631015143),c(638.7388067149211,1664.151603724566,1992.5825553080426,1558.0310637584528),c(1195.0597028859813,994.4320558841919,2493.8953326198753,2665.83936711433),c(1213.0886208174509,1526.1487878059434,1817.6358580512058,2736.537947755733),c(1473.2201509715114,1164.9061226071963,1154.0449374218242,1640.2662328727208),c(1720.473882603094,1712.8584799311386,3158.6927821958557,3960.8953756836304))
targetgene="NON-TARGETING"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1465.4934718580246,1733.1530116838771,3573.738739825869,3842.5713913465543),c(4813.721087702369,6348.129532256637,7768.839887150159,8186.836476282312),c(1859.554106645859,1408.4405036400594,4443.646110323658,3854.6995997411045),c(4167.255601873961,5995.004679758985,12795.239478681073,13677.365159483497),c(798.4235083936513,1175.0533884835654,3354.7537360181727,4179.498936746379),c(3750.0149297456655,3435.8642257386464,6267.314613107858,7996.630671460462))
targetgene="BAX"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(2954.16698105651,2672.7898318356747,8356.419484143811,4172.103687725312),c(3487.3078398871094,2402.8725595242513,8455.354857764918,4965.761812666251),c(1166.7285461365293,2021.3353625727655,3322.1774544600034,1953.5289814051303),c(2019.2388083245892,2419.1081849264424,5467.385921512803,4298.710350965984),c(2093.930039754963,1377.9987060109515,4918.415250810315,2399.018782434223),c(1099.7639938196423,953.8429923787146,3460.928283318874,2018.902982751365))
targetgene="DBR1"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(873.1147398240252,1140.5526845039099,2264.0515682927894,2332.75735120546),c(631.0121276014341,1130.4054186275405,1777.2171383401433,1821.3019289084475),c(1326.4132478152594,777.2805661298887,3304.6827847343197,2202.305158473833),c(2341.1837713865452,2238.4868523270684,5204.362611154249,4227.420150402895),c(2240.736942911215,1418.5877695164288,3869.3383317425933,3268.9958772725763),c(1295.5065313613118,1325.2329234538313,2975.3003823128265,3271.362356959318))
targetgene="TBX3"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(3479.5811607736223,1339.4390956807483,3144.8176993099687,4592.745452043618),c(1416.5578374726072,2613.9356897527327,5004.078806018835,6801.558429655991),c(2730.0932867653887,2410.9903722253466,3872.3546541090905,5941.343063525445),c(700.5522396228166,1020.8149471627521,2002.2347868808336,4324.445817559297),c(909.1725756869643,1313.056204402188,2505.3573576125646,4042.8347348370557),c(2016.6632486200936,1881.3030934788692,3663.625146347485,5467.751316216298))
targetgene="CBLL1"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(3170.5139962341445,5108.133642164308,9159.364498105362,5689.904596809159),c(1632.9048526502418,1725.0351989827818,3457.911960952377,2560.2352110934894),c(6212.250007243507,4596.711441995295,10654.253862941368,6378.845995611785),c(2853.7201525811797,2134.9847403881017,6496.555112961644,4677.938720766313),c(4767.361013021447,2932.559838270729,9385.588675592651,6159.059194705666),c(2073.325562118998,3125.357889921746,8168.200968474386,5035.277153464283))
targetgene="SMAD5"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(1395.953359836642,1869.126374427226,5427.57046627504,4799.812424633502),c(690.250000804834,618.9832184585276,1626.4010200152838,1898.8041386492325),c(2336.032651977554,2622.0535024538285,6683.56709968447,7498.48669740137),c(270.4337689720432,842.2230677386523,775.7981126630768,1392.081675725703),c(1316.111008997277,584.482514478872,2551.8087220566213,1746.4620088152467),c(1952.2742560077024,978.196430482001,4019.551185594153,2998.625573062357))
targetgene="CYCS"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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
targetmat=list(c(826.7546651431035,2565.2288135461604,1770.5812291338495,1819.5270691433914),c(1684.4160467401548,1522.0898814553957,1849.0056106627765,2594.844976512084),c(646.4654858284081,1073.5807297198724,1128.1045650699484,2009.1412540435563),c(2029.5410471425719,1079.669089245694,1489.4599845763116,1859.4614138571546),c(1012.1949638667903,1181.1417480093871,1922.0006119320085,2472.3796527232103),c(937.5037324364164,1014.7265876369305,669.6235653623758,1293.8727687259297),c(2354.0615699090235,1599.2091021158024,1847.195817242878,1631.0961240865972),c(1099.7639938196423,2305.458807111106,1709.0482528573068,1656.8315906799114),c(463.6007468092169,117.70828416588394,939.8860494005239,1077.0440674282374),c(1130.6707102735902,1185.2006543599348,1913.5549093058162,2173.6115922720924),c(1524.7313450614245,744.809315325507,1934.065901397997,2426.8249187534357),c(1462.917912153529,1077.6396360704202,1870.1198672282567,2867.5817604090453),c(1457.7667927445377,836.1347082128307,1050.8867124876203,1594.1198789812609),c(455.87406769573,543.8934509733947,1088.2891098321857,2146.6928858354077),c(252.40485104057365,1714.8879331064124,810.7874521144442,1756.5195474838981),c(1689.567166149146,1317.1151107527357,2092.724457875749,2230.9987246755745),c(303.9160451304866,1475.412458424097,2020.935985553116,952.8038838743072),c(1913.6408604402675,1412.4994099906073,2498.7214484062706,2384.819904313774),c(133.92910463377376,389.65500965258127,679.8790614084662,606.7062296883589),c(800.999068098147,903.1066629968681,946.5219586068177,2085.756033901813),c(579.5009335115211,734.6620494491376,1087.6858453588861,1050.1253609915525),c(880.8414189375121,265.8583659608758,454.25814839447656,955.1703635610487),c(1717.8983228985983,876.723771718308,2020.935985553116,1351.851521051097),c(3149.909518598179,1321.1740171032834,2673.064881189808,2808.1239582796643),c(1764.2583975795199,1883.332546654143,3146.627492729867,3678.101053118018),c(1331.5643672242509,633.1893906854447,2017.919663186619,1365.1629692890183),c(947.8059712543991,1132.4348718028143,2033.001275019105,2378.607895136077),c(592.3787320339994,373.4193842503904,1879.1688343277483,1033.2641932235192),c(1766.8339572840155,2275.0170094819982,2885.4139757912103,2694.8287432769134),c(515.1119408991299,651.4544692629094,585.1665391004545,649.3028640497064),c(839.6324636655818,842.2230677386523,1276.5076255016102,1693.2162158635624),c(798.4235083936513,448.5091517355233,950.7448099199138,669.713751347852),c(692.8255605093298,614.9243121079799,784.2438152892689,1230.5694371055938),c(136.50466433826944,138.00281591862253,576.1175720009629,345.8018442251056),c(229.2248137001128,480.980402539905,739.6022442651106,1313.3962261415475),c(723.7322769632775,635.2188438607185,682.8953837749634,892.4586518623983),c(883.4169786420077,1341.4685488560222,1203.5126242323781,1103.3711539432368),c(1308.3843298837899,685.955173242565,714.8684008598336,1157.2085668166067),c(1828.6473901919112,2155.2792721408405,2893.8596784174024,2781.501061803822),c(1295.5065313613118,549.9818104992163,1007.451670410061,458.5054393061708),c(834.4813442565904,1438.8823012691673,1088.892374305485,1908.861677317884),c(1084.3106355926684,375.4488374256643,1232.4693189507511,1045.9840215397548),c(1725.6250020120851,992.402602708918,2255.6058656665973,2132.789817675801),c(911.74813539146,1674.2988696009352,1037.0116296017334,905.178480178634),c(1223.3908596354336,813.8107232848182,1196.8767150260844,1326.7076743794685),c(1223.3908596354336,1262.3198750203414,2352.731445867807,1228.2029574188523),c(1017.3460832757816,460.68587078716644,1743.434327835375,1675.4676182130008),c(1527.30690476592,618.9832184585276,1012.8810506697558,1454.4975774635109),c(368.30503774287786,1217.6719051643165,1043.647538808027,1293.2811488042444),c(412.0895527193039,229.32820880594628,532.079265450104,1487.036673156207),c(1638.055972059233,1442.9412076197152,2101.7734249752407,1694.6952656677759),c(1614.8759347187722,1213.6129988137689,1882.788421167545,2062.978666916926),c(746.9123143037384,1881.3030934788692,2093.930986822348,1776.930434782044),c(1854.4029872368676,2471.8739674835624,2415.4709510909483,2616.734913614443),c(1568.5158600378506,870.6354121924863,1602.873705556606,2209.996217455743),c(600.1054111474864,874.6943185430341,959.7937770194053,1324.0453847318843),c(2032.1166068470675,1534.266600507039,1593.2214739838148,1848.2206353451325),c(2513.746271587754,1651.9748846729228,2538.536903644034,2553.7273919549502),c(1586.5447779693202,972.1080709561794,1500.922009569001,1300.9722077861543),c(1596.8470167873027,631.1599375101707,1240.311757103644,1214.0040792984032))
targetgene="SAFE"
collabel=c("D0R1","D0R2","D3TOTAL1","D3TOTAL2")

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

