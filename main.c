#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define NUC_NUM 4 //塩基の種数（A, T, G, C）

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列


struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

void make_frequency_matrix(int motif_num);

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

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

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  make_frequency_matrix(seq_num);

  return 0;
}

void make_frequency_matrix(int motif_num)
{
  if(motif_num == 0) return;

  int motif_len = strlen(g_motif[0]);
  int matrix[NUC_NUM][motif_len]; // A=0, T=1, G=2, C=3

  // 頻度表(matrix)初期化
  for(int i = 0; i < NUC_NUM; i++){
    for(int j = 0; j < motif_len; j++){
      matrix[i][j] = 0;
    }
  }

  // 塩基出現頻度カウント
  for(int i = 0; i < motif_num; i++){
    for(int j = 0; j < motif_len; j++){
      char base = g_motif[i][j];
      switch(base){
        case 'A': matrix[0][j]++; break;
        case 'T': matrix[1][j]++; break;
        case 'G': matrix[2][j]++; break;
        case 'C': matrix[3][j]++; break;
      } 
    }
  }

  // 頻度表作成
  printf("A  T  G  C\n");
  for(int j = 0; j < motif_len; j++)
  {
    printf("%d %d %d %d\n", matrix[0][j], matrix[1][j], matrix[2][j], matrix[3][j]);
  }

  float Pi[NUC_NUM][motif_len]; // A=0, T=1, G=2, C=3
  int bunbo=18;

  // 出現確率Pの行列初期化
  for(int i = 0; i < NUC_NUM; i++){
    for(int j = 0; j < motif_len; j++){
      Pi[i][j] = 0;
    }
  }

  // 出現確率Pの行列計算
  for(int i = 0; i < NUC_NUM; i++){
    for(int j = 0; j < motif_len; j++){
      float bunshi=(float)matrix[i][j]+1; // 疑似頻度1を加える
      Pi[i][j]=bunshi/bunbo;
    } 
  }

  //バックグラウンドの出現確率qの計算
  int Q[NUC_NUM];  // 実際の塩基出現頻度
  Q[0]=7519429; Q[1]=7519429; Q[2]=4637676; Q[3]=4637676;
  int totalQ=0;
  for(int i=0;i<NUC_NUM;i++){totalQ+=Q[i];}

  float q[NUC_NUM]; // バックグラウンドの出現確率q
  for(int i = 0; i < NUC_NUM; i++){q[i]=(double)Q[i]/totalQ;}
  
  //オッズスコア行列の計算および出力
  printf("\nオッズスコア\n");\
  double s[NUC_NUM][motif_len];
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < motif_len; j++){
      double x=Pi[i][j]/q[i];
      s[i][j]=log(x);
      printf("%f ", s[i][j]);
    } 
      printf("\n");
  }

//結合部位の探索
int promoter_size=strlen(g_pro[0].seq);
int motif_size=strlen(g_motif[0]);
double hit;
double mean;
double standerd_division, in_brackets, sigma=0;
double total=0, count=0;

for(int l=0; l<2; l++){
for(int i=0; i<MAX_GENE_NUM; i++){
  for(int j=0; j<promoter_size-motif_size+1; j++){
    hit=0;
    for(int k=0; k<motif_size; k++){
      if(g_pro[i].seq[j+k]=='A'){hit += s[0][k];}
      if(g_pro[i].seq[j+k]=='T'){hit += s[1][k];}
      if(g_pro[i].seq[j+k]=='G'){hit += s[2][k];}
      if(g_pro[i].seq[j+k]=='C'){hit += s[3][k];}}

      if(l==0){
           total+=hit; count++;
      if(hit > -1.78){
        printf("\npro: %s\n", g_pro[i].name);
        printf("pos: %d\nhit(", j+1);
        for(int l=0; l<motif_size; l++){
          printf("%c",g_pro[i].seq[j+l]);
        }
        printf(")= %.2f\n\n", hit);
      }}
      if(l==1){
        in_brackets=hit-mean;
        sigma+=in_brackets*in_brackets;
      }   
  }
 }
  if(l==0){mean=total/count; printf("平均値: %f\n", mean);}
  if(l==1){standerd_division=sqrt(sigma/count); printf("標準偏差: %f\n", standerd_division);}
}}