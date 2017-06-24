#include <stdio.h>
#include <string.h>

int main()
{
	char fname[100];
	char Rname[100];
	char line[1024];
	char master[100];
	char method[20];
	char crit[20];
	char acronym[50];
	char home[100];
	char script[100];
	char datapath[100];
	char data[20];
	char src[50];
	char vname[50];
	FILE *f, *h;
	int L, k;
	double v0, lam, eta;

	strcpy(data,"scRNA");
	strcpy(method,"GBC");
	strcpy(crit,"BIC");

	strcpy(master,"/home/cchan40/project/GBC");
	sprintf(home,"%s/GBC",master);
	sprintf(datapath,"%s/dataset/%s/data",master,data);
	sprintf(script,"%s/20170525",home);
	strcpy(src,"DataGBC.R");

	sprintf(fname,"%s%s%s",method,crit,data);
	h = fopen(fname,"w");
	chmod(fname,0755);

	for ( L=6 ; L<=10 ; L++ )
	for ( k=5 ; k<=15 ; k+=5 )
	for ( v0=0.02 ; v0<=0.101 ; v0+=0.02 )
	for ( lam=0.10 ; lam<=0.41 ; lam+=0.05 )
	for ( eta=0.00 ; eta<=0.21 ; eta+=0.1 )
	{
		sprintf(acronym,"%s%s%sL%02dk%02dv0%.2flam%.2feta%.1f",method,crit,data,L,k,v0,lam,eta);
		sprintf(vname,"res%s",acronym);
		sprintf(fname,"%s/%s",script,acronym);

		sprintf(line,"qsub -q fruit.q %s\n",fname);
		fputs(line,h);

		f = fopen(fname,"w");
		sprintf(Rname,"%s.R",fname);
		sprintf(line,"R --vanilla < %s\n",Rname);
		fputs(line,f);
		fclose(f);

		f = fopen(Rname,"w");
		fputs("library(compiler)\n",f);
		fputs("enableJIT(3)\n",f);
		sprintf(line,"source(\"%s/%s\")\n",home,src);
		fputs(line,f);
		sprintf(line,"load(\"%s\")\n",datapath);
		fputs(line,f);

		sprintf(line,"L = %d\n",L);
		fputs(line,f);
		sprintf(line,"k = %d\n",k);
		fputs(line,f);
		sprintf(line,"v0 = %.02f\n",v0);
		fputs(line,f);
		sprintf(line,"lam = %.02f\n",lam);
		fputs(line,f);
		sprintf(line,"eta = %.01f\n",eta);
		fputs(line,f);
		fputs("param = rep(1,p)\n",f);

		sprintf(line,"if ( !file.exist(\"%s/%s\") )\n",script,vname);
		fputs(line,f);
		fputs("{\n",f);
		sprintf(line,"  %s = DataGBC_%s(X,E,L,k,v0,lam,eta,param,intercept=F)\n",vname,crit);
		fputs(line,f);
		sprintf(line,"  save(%s,file=\"%s/%s\")\n",vname,script,vname);
		fputs(line,f);
		fputs("}\n",f);

		fclose(f);

	}
	fclose(h);
}



