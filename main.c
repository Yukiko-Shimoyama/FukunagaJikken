#include<stdio.h>
#include<math.h>
int main(void){
    FILE *fp=fopen("MATa1.txt","r");

    if(fp==NULL){
        printf("cannot open file\n");
        return 0;
    }

    int x, y;
    fscanf(fp, "%d %d\n", &x, &y);
    printf("%d %d\n", x, y);
    fscanf(fp, "%d %d\n", &x, &y);
    printf("%d %d\n", x, y);

    fclose(fp);

    return 0;
}