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
	int K, s, d0, d1, d2, d3;

	strcpy(data,"scRNA");
	strcpy(method,"GBC");
	strcpy(crit,"CCV");

	strcpy(master,"/home/changgee/project/GBC");
	sprintf(home,"%s/GBC",master);
	sprintf(datapath,"%s/dataset/%s/data",master,data);
	sprintf(script,"%s/20170507",home);
	strcpy(src,"DataGBC.R");

	K = 3;
	sprintf(fname,"%s%s%d%s",method,crit,K,data);
	h = fopen(fname,"w");
	chmod(fname,0755);

	sprintf(fname,"%s%s%d%sMERGE",method,crit,K,data);
	m = fopen(fname,"w");
	chmod(fname,0755);

	for ( s=0 ; s<3 ; s++ )
	{
		sprintf(acronym,"%s%s%d%s%d",method,crit,K,data,s);
		sprintf(vname,"res%s",acronym);

		sprintf(fname,"%s/%s",script,acronym);
		sprintf(line,"%s\n",fname);
		fputs(line,h);

		g = fopen(fname,"w");
		chmod(fname,0755);

		for ( d0=1 ; d0<=2 ; d0++ )
		for ( d1=1 ; d1<=3 ; d1++ )
		for ( d2=1 ; d2<=5 ; d2++ )
		for ( d3=1 ; d3<=5 ; d3++ )
			{
				sprintf(fname,"%s/%s%d%d%d%d",script,acronym,d0,d1,d2,d3);
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

				sprintf(line,"%s = DataGBC_%s(X,E,L,k,v0,lam,eta,param,intercept=F,smoothing=smoothing,fold=%d,run=c(%d,%d,%d,%d))\n",vname,crit,K,d0,d1,d2,d3);
				fputs(line,f);

				sprintf(line,"save(%s,file=\"%s/%s%d%d%d%d\")\n",vname,script,vname,d0,d1,d2,d3);
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

		sprintf(line,"for ( k in 1:2 )\n");
		fputs(line,g);
		sprintf(line,"for ( l in 1:3 )\n");
		fputs(line,g);
		sprintf(line,"for ( i in 1:5 )\n");
		fputs(line,g);
		sprintf(line,"for ( j in 1:5 )\n");
		fputs(line,g);
		fputs("{\n",g);
		sprintf(line,"  fname = sprintf(\"%s/%s%%d%%d%%d%%d\",k,l,i,j)\n",script,vname);
		fputs(line,g);
		fputs("  load(fname)\n",g);
		fputs("  if ( k==1 & l==1 & i==1 & j==1 )\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp = %s\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("  else\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp$CCV = tmp$CCV + %s$CCV\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("}\n",g);
		fputs("idx = which.min(apply(tmp$CCV,c(1,2),sum))\n",g);
		fputs("i = (idx-1)%%5 + 1\n",g);
		fputs("j = (idx-1)%/%5 + 1\n",g);
		fputs("tmp$opt_v0 = tmp$v0[i]\n",g);
		fputs("tmp$opt_lam = tmp$lam[j]\n",g);
		sprintf(line,"fname = sprintf(\"%s/%s11%%d%%d\",i,j)\n",script,vname);
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



