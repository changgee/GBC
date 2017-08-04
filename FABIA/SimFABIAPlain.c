#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define HPC 0
#define LPC 1
#define Emory 2

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
	int s, batch_size, batch, R, where;

	R = 100;
	batch_size = 10;

	strcpy(method,"FABIA");
	strcpy(crit,"Plain");

	if ( access("/home/cchan40",X_OK) == 0 )
	{
		where = Emory;
		strcpy(master,"/home/cchan40/project/GBC");
	}
	else if ( access("/project/qlonglab/changgee",X_OK) == 0 )
	{
		where = HPC;
		strcpy(master,"/home/changgee/project/GBC");
	}
	else
	{
		where = LPC;
		strcpy(master,"/home/changgee/project/GBC");
	}
	sprintf(home,"%s/FABIA",master);
	sprintf(script,"%s/Sim%s",home,crit);
	strcpy(src,"SimFABIA.R");

	sprintf(fname,"%s%s",method,crit);
	h = fopen(fname,"w");
	chmod(fname,0755);

	sprintf(fname,"%s%sMERGE",method,crit);
	m = fopen(fname,"w");
	chmod(fname,0755);

	for ( s=0 ; s<16 ; s++ )
	{
		sprintf(acronym,"%s%s%02d",method,crit,s+1);
		sprintf(vname,"res%s",acronym);

		sprintf(fname,"%s/%s",script,acronym);
		sprintf(line,"%s\n",fname);
		fputs(line,h);

		g = fopen(fname,"w");
		chmod(fname,0755);

		for ( batch=0 ; batch<R ; batch+=batch_size )
		{
			sprintf(fname,"%s/%s%03d",script,acronym,batch+1);
			if ( where == Emory )
				sprintf(line,"qsub -q fruit.q %s\n",fname);
			else if ( where == HPC )
				sprintf(line,"bsub -q qlonglab -e %s.e -o %s.o < %s\n",fname,fname,fname);
			else
				sprintf(line,"bsub -q cceb_normal -e %s.e -o %s.o < %s\n",fname,fname,fname);
			fputs(line,g);

			f = fopen(fname,"w");
			if ( where == LPC )
				fputs("module load R\n",f);
			else if ( where == HPC )
			{
				fputs("source /etc/profile.d/modules.sh\n",f);
				fputs("module load R-3.3.1\n",f);
			}
			sprintf(Rname,"%s.R",fname);
			sprintf(line,"R --vanilla < %s\n",Rname);
			fputs(line,f);
			fclose(f);

			f = fopen(Rname,"w");
			fputs("library(compiler)\n",f);
			fputs("enableJIT(3)\n",f);
			sprintf(line,"source(\"%s/%s\")\n",home,src);
			fputs(line,f);

			if ( (s/8)%2 == 0 )
			{
				fputs("type = 0\n",f);
				fputs("param = 4\n",f);
			}
			else
			{
				fputs("type = NULL\n",f);
				fputs("param = NULL\n",f);
			}

			if ( (s/2)%4 == 0 )
				fputs("L = 4\n",f);
			else if ( (s/2)%4 == 1 )
				fputs("L = 5\n",f);
			else if ( (s/2)%4 == 2 )
				fputs("L = 10\n",f);
			else
				fputs("L = 20\n",f);

			if ( s%2 == 0 )
				fputs("overlap = 0\n",f);
			else
				fputs("overlap = 15\n",f);

			fputs("p = 1000\n",f);
			fputs("n = 300\n",f);

			fputs("seed = 100\n",f);
			fputs("thrW = 2:6/2\n",f);
			fputs("thrZ = 1:5/10\n",f);

			sprintf(line,"%s = SimFABIA_%s(%d,seed,p,n,type,param,overlap,L,thrW,thrZ,batch=%d)\n",vname,crit,batch_size,batch);
			fputs(line,f);
			sprintf(line,"save(%s,file=\"%s/%s%03d\")\n",vname,script,vname,batch+1);
			fputs(line,f);
			fclose(f);
		}

		fclose(g);


		sprintf(fname,"%s/%sMERGE",script,acronym);
		if ( where == Emory )
			sprintf(line,"qsub -q fruit.q %s\n",fname);
		else if ( where == HPC )
			sprintf(line,"bsub -q qlonglab -e %s.e -o %s.o < %s\n",fname,fname,fname);
		else
			sprintf(line,"bsub -q cceb_normal -e %s.e -o %s.o < %s\n",fname,fname,fname);
		fputs(line,m);

		f = fopen(fname,"w");
		if ( where == LPC )
			fputs("module load R\n",f);
		else if ( where == HPC )
		{
			fputs("source /etc/profile.d/modules.sh\n",f);
			fputs("module load R-3.3.1\n",f);
		}
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
//		sprintf(line,"    tmp$fits = c(tmp$fits,%s$fits)\n",vname);
//		fputs(line,g);
		sprintf(line,"    tmp$CE = abind(tmp$CE,%s$CE)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$FP = abind(tmp$FP,%s$FP)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$FN = abind(tmp$FN,%s$FN)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$SEN = abind(tmp$SEN,%s$SEN)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$SPE = abind(tmp$SPE,%s$SPE)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$MCC = abind(tmp$MCC,%s$MCC)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$CS = abind(tmp$CS,%s$CS)\n",vname);
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



