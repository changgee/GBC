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
	char acronym[20];
	char home[100];
	char script[100];
	char datapath[100];
	char data[20];
	char src[50];
	char vname[50];
	FILE *f, *g, *h, *m;
	int s, d1, d2;

	strcpy(data,"scRNA");
	strcpy(method,"GBC");
	strcpy(crit,"BIC");

	strcpy(master,"/home/changgee/project/GBC");
	sprintf(home,"%s/GBC",master);
	sprintf(datapath,"%s/dataset/%s/data",master,data);
	sprintf(script,"%s/20170507",home);
	strcpy(src,"DataGBC.R");

	sprintf(fname,"%s%s%s",method,crit,data);
	h = fopen(fname,"w");
	chmod(fname,0755);

	sprintf(fname,"%s%s%sMERGE",method,crit,data);
	m = fopen(fname,"w");
	chmod(fname,0755);

	for ( s=0 ; s<3 ; s++ )
	{
		sprintf(acronym,"%s%s%s%d",method,crit,data,s);
		sprintf(vname,"res%s",acronym);

		sprintf(fname,"%s/%s",script,acronym);
		sprintf(line,"%s\n",fname);
		fputs(line,h);

		g = fopen(fname,"w");
		chmod(fname,0755);

		for ( d1=1 ; d1<=5 ; d1++ )
			for ( d2=1 ; d2<=5 ; d2++ )
			{
				sprintf(fname,"%s/%s%d%d",script,acronym,d1,d2);
				sprintf(line,"bsub < %s\n",fname);
				fputs(line,g);

				f = fopen(fname,"w");
				fputs("module load R\n",f);
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
				fputs("set.seed(100)\n",f);

				fputs("L = 9\n",f);
				fputs("k = 10\n",f);
				fputs("param = rep(1,p)\n",f);

				fputs("v0 = (4:8/20)^2\n",f);
				fputs("lam = 6:10/20\n",f);

				if ( s == 0 )
				{
					fputs("eta = 0\n",f);
					fputs("smoothing = \"Ising\"\n",f);
				}
				else if ( s == 1 )
				{
					fputs("eta = 0.2\n",f);
					fputs("smoothing = \"Ising\"\n",f);
				}
				else if ( s == 2 )
				{
					fputs("eta = 0.2\n",f);
					fputs("smoothing = \"MRF\"\n",f);
				}

				sprintf(line,"%s = DataGBC_%s(X,E,L,k,v0,lam,eta,param,intercept=F,smoothing=smoothing,run=c(%d,%d))\n",vname,crit,d1,d2);
				fputs(line,f);

				sprintf(line,"save(%s,file=\"%s/%s%d%d\")\n",vname,script,vname,d1,d2);
				fputs(line,f);
				fclose(f);

			}

		fclose(g);


		sprintf(fname,"%s/%sMERGE",script,acronym);
		sprintf(line,"bsub < %s\n",fname);
		fputs(line,m);

		f = fopen(fname,"w");
		fputs("module load R\n",f);
		sprintf(Rname,"%s.R",fname);
		sprintf(line,"R --vanilla < %s\n",Rname);
		fputs(line,f);
		fclose(f);

		g = fopen(Rname,"w");

		sprintf(line,"for ( i in 1:5 )\n");
		fputs(line,g);
		sprintf(line,"  for ( j in 1:5 )\n");
		fputs(line,g);
		fputs("  {\n",g);
		sprintf(line,"    fname = sprintf(\"%s/%s%%d%%d\",i,j)\n",script,vname);
		fputs(line,g);
		fputs("    load(fname)\n",g);
		fputs("    if ( i==1 & j==1 )\n",g);
		fputs("    {\n",g);
		sprintf(line,"      tmp = %s\n",vname);
		fputs(line,g);
		fputs("    }\n",g);
		fputs("    else\n",g);
		fputs("    {\n",g);
		sprintf(line,"      tmp$BIC = tmp$BIC + %s$BIC\n",vname);
		fputs(line,g);
		fputs("    }\n",g);
		fputs("  }\n",g);
		fputs("idx = which.min(tmp$BIC)\n",g);
		fputs("i = (idx-1)%%5 + 1\n",g);
		fputs("j = (idx-1)%/%5 + 1\n",g);
		fputs("tmp$opt_v0 = tmp$v0[i]\n",g);
		fputs("tmp$opt_lam = tmp$lam[j]\n",g);
		sprintf(line,"fname = sprintf(\"%s/%s%%d%%d\",i,j)\n",script,vname);
		fputs(line,g);
		fputs("load(fname)\n",g);
		sprintf(line,"tmp$opt_fit = %s$opt_fit\n",vname);
		fputs(line,g);
		sprintf(line,"tmp$opt_biclus = %s$opt_biclus\n",vname);
		fputs(line,g);
		sprintf(line,"tmp$time = %s$time\n",vname);
		fputs(line,g);
		sprintf(line,"%s = tmp\n",vname);
		fputs(line,g);
		sprintf(line,"save(%s,file=\"%s/%s\")\n",vname,home,vname);
		fputs(line,g);

		fclose(g);

	}
	fclose(h);
	fclose(m);
}



