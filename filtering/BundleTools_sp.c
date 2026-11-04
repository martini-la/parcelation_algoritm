#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "BundleTools_sp.h"

struct bundle read_bundle(char* bunfile)
{
    int i, j;
    struct bundle fas;

    FILE *fb;
    fb = fopen(bunfile, "r");
    if (fb == NULL) {fputs ("File error",stderr); exit (1);}

    char buffer[100];   //buffer para la linea
    fgets(buffer,100,fb);//Lecturas para llegar la liena 5 del archivo
    fgets(buffer,100,fb);
    fgets(buffer,100,fb);
    fgets(buffer,100,fb);
    fgets(buffer,100,fb);
    int index_i=0;//Indice posición inicial de los numeros
    for(i=0;i<sizeof(buffer);i++)
    {
        if(buffer[i]>=48 && buffer[i]<=57)//Detecta cuando hay un número en la línea leída y guarda el índice
        {
            index_i=i;
            break;
        }
    }
    int index_f=index_i;//Indice posición final de los números
    int32_t nfibers=0;   // variable para numero de fibras
    while(buffer[index_f]!=44)//Mientras no se detecte una coma aumenta el otro indice de la pos final
    {
        index_f++;
    }
    for(i=0;i<(index_f-index_i);i++)
    {
        nfibers+=(buffer[index_i+i]-48)*pow(10,(index_f-(index_i+i))-1);
    }
    fclose(fb);

    fas.nfibers = nfibers;

    char* bunfileb;
    bunfileb = masdata(bunfile);

    FILE *fp;
    fp = fopen(bunfileb, "rb");
    if (fp == NULL) {fputs ("File error",stderr); exit (1);}

    float** points; // puntero a cada fibra
    points = (float**) malloc (nfibers*sizeof(float*));

    int32_t* npoints;//Puntero para la cantidad de puntos de cada fibra
    npoints=(int32_t*) malloc(nfibers*sizeof(int32_t));//Memoria para cada fibra

    for ( i = 0; i < nfibers; i++ )//Itera en cantidad de fibras
    {
       fread(npoints+i,sizeof(int32_t),1,fp);//Lee primer elemento(cantidad de funtos)

       points[i]=(float*) malloc((*(npoints+i))*3*sizeof(float));//Asigna memoria para toda una fibra
       if (points[i] == NULL) {fputs ("Memory error",stderr); exit (2);}
       fread( points[i], sizeof(float), *(npoints+i)*3, fp );//Lee todos los puntos de la fibra

    }
    fas.npoints=npoints;
    fas.points=points;
    fclose(fp);
    return fas;
}

void write_bundle(char* outfile, int32_t nfibers, int32_t* npoints, float** points)
{
    char par1[] = "attributes = {\n    'binary' : 1,\n    'bundles' : ['fibers', 0],\n    'byte_order' : 'DCBA',\n    'curves_count' : ";
    char par2[] = ",\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }";

    int32_t len=strlen(int2string(nfibers))+1;

    FILE *fw;
    fw = fopen(outfile,"w");

    fwrite(par1,sizeof(char),strlen(par1),fw);
    fwrite(int2string(nfibers),sizeof(char),len,fw);
    fseek(fw, -1 , SEEK_CUR);
    fwrite(par2,sizeof(char),strlen(par2),fw);

    fclose(fw);
    char* outfileb;
    outfileb=masdata(outfile);

    FILE *fwb;
    fwb = fopen(outfileb,"wb");
    int i;
    for(i=0;i<nfibers;i++)
    {
        fwrite(npoints+i,sizeof(int32_t),1,fwb);
        fwrite(points[i],sizeof(float),*(npoints+i)*3,fwb);
    }

    fclose(fwb);
}

char* masdata(char* bunfile)
{
    int cont=0,i;
    char data[5]="data";
    char* bunfileb=malloc(sizeof(char)*(strlen(bunfile)+strlen(data)));

    while(bunfile[cont]!=NULL)
    {
        bunfileb[cont]=bunfile[cont];//Iguala string excepto el nulo
        cont++;
    }
    for(i=0;i<sizeof(data);i++)
    {
        bunfileb[cont]=data[i];//Le agrega data al final
        cont++;
    }
    return bunfileb;
}



char* int2string(int32_t si)
{
    int32_t exp,base,i;
    for(exp=0;;exp++)
    {
        base=pow(10,exp);
        if(base==99){base++;}
        if(si<base)
        {
            exp--;
            base=pow(10,exp);
            if(base==99){base++;}
            break;
        }
    }

    char* cad=(char*)malloc(sizeof(char)*(exp+1));
    for(i=0;i<exp+1;i++)
    {
        cad[i]=(char)((si/base)%10 +48);
        base=pow(10,exp-(i+1));
        if(base==99){base++;}
    }
    cad[i]=NULL;
    return cad;
}
