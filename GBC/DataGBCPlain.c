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
	char acronym[50];
	char home[100];
	char script[100];
	char datapath[100];
	char data[20];
	char ver[20];
	char src[50];
	char vname[50];
	FILE *f, *h;
	int L, k, where;
	double eta;

	strcpy(data,"NCI60");
	strcpy(ver,"1000");
	strcpy(method,"GBC");
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
	sprintf(home,"%s/GBC",master);
	sprintf(datapath,"%s/dataset/%s/data%s",master,data,ver);
	sprintf(script,"%s/%s%s",home,data,crit);
	strcpy(src,"DataGBC.R");

	sprintf(fname,"%s%s%s%s",method,crit,data,ver);
	h = fopen(fname,"w");
	chmod(fname,0755);

	for ( L=10 ; L<=10 ; L++ )
	for ( k=10 ; k<=20 ; k+=5 )
	for ( eta=0.00 ; eta<=0.21 ; eta+=0.05 )
	{
		sprintf(acronym,"%s%s_%s%s_L%02d_k%02d_eta%.2f",method,crit,data,ver,L,k,eta);
		sprintf(vname,"res%s",acronym);
		sprintf(fname,"%s/%s",script,acronym);
		if ( where == Emory )
			sprintf(line,"qsub -q fruit.q %s\n",fname);
		else if ( where == HPC )
			sprintf(line,"bsub -q qlonglab -e %s.e -o %s.o < %s\n",fname,fname,fname);
		else
			sprintf(line,"bsub -q cceb_normal -e %s.e -o %s.o < %s\n",fname,fname,fname);
		fputs(line,h);

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
		sprintf(line,"load(\"%s\")\n",datapath);
		fputs(line,f);

		sprintf(line,"L = %d\n",L);
		fputs(line,f);
		sprintf(line,"k = %d\n",k);
		fputs(line,f);
		sprintf(line,"v0 = 1:9/1000\n");
		fputs(line,f);
		sprintf(line,"lam = 1:10/1000\n");
		fputs(line,f);
		sprintf(line,"eta = %.02f\n",eta);
		fputs(line,f);

		sprintf(line,"dpath = '%s'\n",datapath);
		fputs(line,f);
		sprintf(line,"opath = '%s'\n",script);
		fputs(line,f);
		sprintf(line,"name = '%s%s'\n",data,ver);
		fputs(line,f);

		sprintf(line,"DataGBC_%s(dpath,opath,name,L,k,v0,lam,eta)\n",crit);
		fputs(line,f);

		fclose(f);

	}
	fclose(h);
}





