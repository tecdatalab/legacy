runReg	<-	function(data)	
{
	data$hit	<-	(data$rmsd<=2.5)*1;
	mod		<-	glm(hit~v1+v2+v3+v4+v5+v6+v7+v8+v9+v10+v11+v12+id,family="binomial",data);
	smod	<-	step(mod,trace=F);
	kk		<-	rep(0,12);

	cf		<-	coefficients(smod);
	cf		<-	cf[names(cf) %in% vnames];
	kk[vnames[2:13] %in% names(cf)]	<-	-cf*1000;

	print(kk);

}

vnames	<-	c('id',paste("v",1:12,sep=""),'rmsd');

out	<-	NULL;

data	<-	read.table('/fenix/corgi2d7/yang41/sigscore_cv/Regression/soro_zdockbm2_pdbh_smalldat.dat',sep=",",header=F);
names(data)	<-	vnames;
out	<-	rbind(out,runReg(data));

data	<-	subset(data,rmsd<=9);
cf		<-	runReg(data);
out	<-	rbind(out,runReg(data));

data	<-	read.table('/fenix/corgi2d7/yang41/sigscore_cv/Regression/soro_zdbm2_c2h_sm.dat',sep=",",header=F);
names(data)	<-	vnames;
out	<-	rbind(out,runReg(data));

row.names(out)	<-	c('sm','sm_c2h','c2h_sm');
write.table(as.data.frame(out),file='logit12.weights',sep=" ",row.names=T,col.names=F,quote=F);
