#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h> 

#define TASK_SIZE 1000

void writeToFile(float* array, int lines, int type){ 
    //quick 
    if(type == 0) {
        FILE* file = fopen("outputquick.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    FILE* file = fopen("outputselect.txt", "w");
    int i;
    for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);

}

int countLines(char* filePath){ 
    FILE* file = fopen(filePath, "r");
    char c;
    int count = 0;
    for (c = getc(file); c != EOF; c = getc(file))
            if (c == '\n') count++;

    return count;
}

float* readFromFile (char* filePath, int lines) { 
    float* array = malloc(lines*(sizeof(float)));
    int i = 0;
    float num;
    FILE* file = fopen(filePath, "r");
    while(fscanf(file, "%f", &num) > 0){
        array[i] = num;
        i++;
    }
    fclose(file);
    return array;
}

void swap(float *a, float *b) {
    float t = *a;
    *a = *b;
    *b = t;
}

// ------ FUNÇÕES IMPLEMENTAÇÃO QUICK SORT ------------- //

int partition(float array[], int low, int high) {
    // pivot = maior elemento do array
    float pivot = array[high];
    // index do menor elemento 
    // i = indica elemento a esquerda do pivot 
    int i = (low - 1);
    //compara todos elementos do array com o pivot
    for (int j = low; j < high; j++) {
        // Se o elemento for menor que o pivot
        if (array[j] <= pivot) {
            i++;
            swap(&array[i], &array[j]);
        }
    }
    // Inverte o pivot com o maior elemento
    swap(&array[i + 1], &array[high]);
    // Retorna posicao do pivot //TODO: ver se esse comentario ta certo
    return (i + 1);
}

int quickSort(float array[], int low, int high){ 
    if (low < high) {
        //faz primeira iteracao, encontrando onde particionar
        int pi = partition(array, low, high);
        //chamada recursiva para array da esquerda 
        quickSort(array, low, pi - 1);
        //chamada recursiva para array da direita
        quickSort(array, pi + 1, high);
    }
}

// --------------------------------------------------------- // 


// ------- FUNÇÕES PARA IMPLEMENTAÇÃO SELECTION SORT ------- //

void selectionSort(float array[], int n)
{
    int i, j, min;
 
    //iteracao na parte nao ordenada do array
    for (i = 0; i < n-1; i++)
    {
        //Encontra o menor valor no array 
        min = i;
        for (j = i+1; j < n; j++) {
            if (array[j] < array[min])
                min = j;
        }
 
        //Inverte menor valor com o primeiro elemento
        if(min != i) swap(&array[min], &array[i]);
    }
}

// ------------------------------------------------------- // 



// ------- FUNÇÕES PARA IMPLEMENTAÇÃO QUICK SORT PARALELO ------- //


int parallel_partion(float * array, int low, int high)
{
    float* lt = malloc (sizeof(float) * (high - low));
    float* gt = malloc (sizeof(float) * (high - low));
    int i;
    int j;
    float key = array[high];
    int lt_n = 0;
    int gt_n = 0;

    for(i = low; i < high; i++){
        if(array[i] < array[high]) lt[lt_n++] = array[i];
        else gt[gt_n++] = array[i];
    }   

    for(i = 0; i < lt_n; i++) array[low + i] = lt[i];
    array[low + lt_n] = key;
    for(j = 0; j < gt_n; j++) array[low + lt_n + j + 1] = gt[j];

    return low + lt_n;
}

void parallel_quickSort(float * array, int low, int high)
{
    int p;

    if(low < high){ 
        p = parallel_partion(array, low, high); 
        #pragma omp task shared(array) if(high - p > TASK_SIZE) 
        parallel_quickSort(array, low, p - 1); 
        #pragma omp task shared(array) if(high - p > TASK_SIZE)
        parallel_quickSort(array, p + 1, high); 
    }
}

// ----------------------------------------------------------------------------------- // 


int main (int argc, char** argv){ 

    printf("Started app");
    int lines = countLines(argv[1]);
    printf("Finished couting lines");
    double start, end;
/*
    //QUICK SORT
    float* arrayQuick = readFromFile(argv[1], lines);
    //------------------------------------------------------------
    double time_quick;
    start = clock();
    quickSort(arrayQuick, 0, lines - 1);
    end = clock();
    time_quick = ((double) (end - start)) / CLOCKS_PER_SEC; 
    //-------------------------------------------------------------
    writeToFile (arrayQuick, lines, 0);


    //SELECTION SORT
    float* arraySelect = readFromFile(argv[1], lines);
    //------------------------------------------------------------
    double time_select;
    start = clock();
    selectionSort(arraySelect, lines);
    end = clock();
    time_select = ((double) (end - start)) / CLOCKS_PER_SEC; 
    //------------------------------------------------------------
    writeToFile (arraySelect, lines, 1);




    //QUICK SORT
    float* arrayQuick = readFromFile(argv[1], lines);
    //------------------------------------------------------------
    start = omp_get_wtime();
    quickSort(arrayQuick, 0, lines - 1);
    end = omp_get_wtime();
    printf("QuickSort: %lf\n", (end - start));
    //-------------------------------------------------------------
*/
    //PARALLEL QUICK SORT


    omp_set_dynamic(0);
    omp_set_num_threads(1);
    printf("Started read from file");
    float* arrayPQuick = readFromFile(argv[1], lines);
    printf("finished read from file and started parallel");
    //------------------------------------------------------------
    start = omp_get_wtime();
    #pragma omp parallel 
    {
        #pragma omp single
        parallel_quickSort(arrayPQuick, 0, lines -1);
    }
    end = omp_get_wtime();
    printf("Parallel QuickSort: %lf", (end - start));
    //------------------------------------------------------------

    return 0;
}


