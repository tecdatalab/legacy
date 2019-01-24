runReg	<-	function(dat)	
{
	dat$hit	<-	(dat$rmsd<=2.5)*1;
	mod		<-	lm(rmsd~v2+v3+v5+v6+v7+v8+v9+v11+v12+id,dat);
	smod	<-	step(mod,trace=F);
	kk		<-	rep(0,12);

	cf		<-	coefficients(smod);
	cf		<-	cf[names(cf) %in% vnames];
	kk[vnames[2:13] %in% names(cf)]	<-	cf*1000;
	print(kk);
}

vnames	<-	c('id',paste("v",1:12,sep=""),'rmsd');

out			<-	NULL;
data		<-	read.table('soro_zdbm2_sm_500.dat',sep=",",header=F);
names(data)	<-	vnames;
out			<-	rbind(out,runReg(data));

data		<-	read.table('soro_zdbm2_c2hsm_500.dat',sep=",",header=F);
names(data)	<-	vnames;
out			<-	rbind(out,runReg(data));

row.names(out)	<-	c('lin_sm','lin_c2hsm');
write.table(as.data.frame(out),file='lin_n9.weights',sep=" ",row.names=T,col.names=F,quote=F);
