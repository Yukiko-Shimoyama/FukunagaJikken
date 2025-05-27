#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>

#define A_frequency 7519429
#define T_frequency 7519429
#define G_frequency 4637676
#define C_frequency 4637676
#define total_frequency 24314210 //バックグラウンド出現確率の合計
#define BUFSIZE 1024  //ファイルから読み込む一行の最大文字数
#define MAX_GENE_NUM 8 //与えられるプロモータ領域の最大遺伝子数

// ランダムに塩基を選ぶ
char choose_base(){
    int x=rand() %total_frequency;
    if(x<A_frequency){return 'A';}
    else if(x<A_frequency+T_frequency){return 'T';}
    else if(x<A_frequency+T_frequency+G_frequency){return 'G';}
    else if(x<A_frequency+T_frequency+G_frequency+C_frequency){return 'C';}
}

//遺伝子のプロモータ領域を保存する構造体
struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; 

//プロモーター配列をファイルから読み込む
int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

//ランダム配列生成、ファイルへの出力
void generate_random_seq(int gene_num, char *result_file){
    FILE *result=fopen(result_file, "w");
    if (result == NULL) {
        perror("output error\n");
        exit(1);
    }
    
    for(int i=0; i<gene_num; i++){
      int length=strlen(g_pro[i].seq);
      for(int j=0; j<length; j++){
          char base=choose_base(); //choose_base関数で選んだ塩基を入力
          fputc(base, result);
        }
    fputc('\n',result);
    }
    fclose(result);
}

int main(int argc, char *argv[]){
  srand((unsigned int)time(NULL));

  char *promoter_file=argv[1];
  char *result_file="random_seq.txt";  

  int gene_num=read_promoter(promoter_file);
  generate_random_seq(gene_num, result_file);

  return 0;

}