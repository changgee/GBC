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
	char src[50];
	char vname[50];
	FILE *f, *g, *h, *m;
	int K, s, batch_size, batch, R;

	K = 3;
	R = 100;
	batch_size = 1;

	strcpy(method,"GBC");
	strcpy(crit,"CCV");

	strcpy(master,"/home/cchan40/project/GBC");
	sprintf(home,"%s/GBC",master);
	sprintf(script,"%s/20170510",home);
	strcpy(src,"SimGBC.R");

	sprintf(fname,"%s%s%d",method,crit,K);
	h = fopen(fname,"w");
	chmod(fname,0755);

	sprintf(fname,"%s%s%dMERGE",method,crit,K);
	m = fopen(fname,"w");
	chmod(fname,0755);

	for ( s=0 ; s<2 ; s++ )
	{
		sprintf(acronym,"%s%s%d%02d",method,crit,K,s+1);
		sprintf(vname,"res%s",acronym);

		sprintf(fname,"%s/%s",script,acronym);
		sprintf(line,"%s\n",fname);
		fputs(line,h);

		g = fopen(fname,"w");
		chmod(fname,0755);

		for ( batch=0 ; batch<R ; batch+=batch_size )
		{
			sprintf(fname,"%s/%s%03d",script,acronym,batch+1);
			sprintf(line,"qsub -q fruit.q %s\n",fname);
			fputs(line,g);

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

			if ( (s/6)%3 == 0 )
			{
				fputs("eta = 0\n",f);
				fputs("smoothing = \"Ising\"\n",f);
			}
			else if ( (s/6)%3 == 1 )
			{
				fputs("eta = 0.1\n",f);
				fputs("smoothing = \"Ising\"\n",f);
			}
			else
			{
				fputs("eta = 0.1\n",f);
				fputs("smoothing = \"MRF\"\n",f);
			}

			if ( (s/2)%3 <= 1 )
				fputs("L = 4\n",f);
			else
				fputs("L = 5\n",f);

			if ( (s/2)%3 == 0 )
			{
				fputs("sigma2 = 9\n",f);
				fputs("seed = 100\n",f);
			}
			else
			{
				fputs("sigma2 = 25\n",f);
				fputs("seed = 200\n",f);
			}

			if ( s%2 == 0 )
				fputs("overlap = 0\n",f);
			else
				fputs("overlap = 15\n",f);

			fputs("k = 10\n",f);
			fputs("v0 = 3:7/30\n",f);
			fputs("lam = 5:9/8\n",f);

			sprintf(line,"%s = SimGBC_%s(%d,seed,overlap,sigma2,L,k,v0,lam,eta,intercept=F,smoothing=smoothing,fold=%d,batch=%d)\n",vname,crit,batch_size,K,batch);
			fputs(line,f);

			sprintf(line,"save(%s,file=\"%s/%s%03d\")\n",vname,script,vname,batch+1);
			fputs(line,f);
			fclose(f);
		}

		fclose(g);


		sprintf(fname,"%s/%sMERGE",script,acronym);
		sprintf(line,"qsub -q fruit.q %s\n",fname);
		fputs(line,m);

		f = fopen(fname,"w");
		sprintf(Rname,"%s.R",fname);
		sprintf(line,"R --vanilla < %s\n",Rname);
		fputs(line,f);
		fclose(f);

		g = fopen(Rname,"w");

		fputs("library(abind)\n",g);
		sprintf(line,"for ( i in 1:%d )\n",R/batch_size);
		fputs(line,g);
		fputs("{\n",g);
		sprintf(line,"  fname = sprintf(\"%s/%s%%03d\",(i-1)*%d+1)\n",script,vname,batch_size);
		fputs(line,g);
		fputs("  load(fname)\n",g);
		fputs("  if ( i == 1 )\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp = %s\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("  else\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp$S = c(tmp$S,%s$S)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$Shat = c(tmp$Shat,%s$Shat)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$CE = c(tmp$CE,%s$CE)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$FP = c(tmp$FP,%s$FP)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$FN = c(tmp$FN,%s$FN)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$SEN = c(tmp$SEN,%s$SEN)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$SPE = c(tmp$SPE,%s$SPE)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$MCC = c(tmp$MCC,%s$MCC)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$CCV = abind(tmp$CCV,%s$CCV)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$opt_v0 = c(tmp$opt_v0,%s$opt_v0)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$opt_lam = c(tmp$opt_lam,%s$opt_lam)\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("}\n",g);
		sprintf(line,"%s = tmp\n",vname);
		fputs(line,g);
		sprintf(line,"save(%s,file=\"%s/%s\")\n",vname,home,vname);
		fputs(line,g);

		fclose(g);

	}
	fclose(h);
	fclose(m);
}



