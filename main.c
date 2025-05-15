#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

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
  void make_frequency_matrix(seq_num);

  return 0;
}


/*void f(char *filename){
    for(int i=0; i<seq_num; i++){
        for(j=0; j<motif_len; j++){
            char bace=g_motif[i][j];
             switch(bace[i]){
                case 'A' : matrix[A][i]++; break;
                case 'C' : matrix[C][i]++; break;
                case 'G' : matrix[G][i]++; break;
                case 'T' : matrix[Y][i]++; break;
                default : fprint(stderr, "unknown char %c/n", bace);
            }
    }
    }
    
    for (i=0; i<4; k++){

        for(int j=0; j<MAX_GENE_NUM; j++){
            printf(matrix[j][i]);}
        printf("\n");

        }

    }
  }   
*/


void make_frequency_matrix(int motif_num)
{
  if(motif_num == 0) return;

  int motif_len = strlen(g_motif[0]);
  int freq[motif_len][4]; // A=0, T=1, G=2, C=3

  // 初期化
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < motif_len; j++){
      freq[j][i] = 0;
    }
  }

  // カウント
  enum dna {A, C, G, T};

  for(int i = 0; i < motif_num; i++){
    for(int j = 0; j < motif_len; j++){
      char base = g_motif[i][j];
      switch(base){
        case 'A': freq[A][j]++; break;
        case 'T': freq[C][j]++; break;
        case 'G': freq[G][j]++; break;
        case 'C': freq[T][j]++; break;
        default:
          printf("Unexpected base '%c' at motif %d position %d\n", base, i, j);
          break;
      }
    }
  }

  // 出力
  printf("Position\tA\tT\tG\tC\n");
  for(int j = 0; j < motif_len; j++)
  {
    printf("%d\t\t%d\t%d\t%d\t%d\n", j + 1, freq[0][j], freq[1][j], freq[2][j], freq[3][j]);
  }

}