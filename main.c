#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h> 

#define TASK_SIZE 1000
struct Compare { float val; int index; };
#pragma omp declare reduction(maximum : struct Compare : omp_out = omp_in.val > omp_out.val ? omp_in : omp_out)

void writeToFile(float* array, int lines, char* type){ 
    if(!strcmp(type, "Q")) {
        FILE* file = fopen("outputq.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    if(!strcmp(type, "PQ")) {
        FILE* file = fopen("outputpq.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    if(!strcmp(type, "S")) {
        FILE* file = fopen("outputs.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    if(!strcmp(type, "PS")) {
        FILE* file = fopen("outputps.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }
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

// -- FUNÇÕES PARA IMPLEMENTAÇÃO SELECTION SORT PARALELO -- //


void parallel_selectionSort(float* array, int lines){
	
    int startpos;
	for(startpos = 0; startpos < lines; startpos++){
		struct Compare max;
        max.val = array[startpos];
        max.index = startpos;

        #pragma omp parallel for reduction(maximum:max)
		for(int i=startpos +1; i< lines; ++i){
			if(array[i] > max.val){
				max.val = array[i];
				max.index = i;
			}
		}

		swap(&array[startpos], &array[max.index]);
	}
}


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

    double start, end;
    int lines = countLines(argv[1]);
    float* array = readFromFile(argv[1], lines);

    omp_set_dynamic(0);
    omp_set_num_threads(4);

    if(!strcmp(argv[2], "Q")) {
        printf("Execução com QuickSort\n");
        start = omp_get_wtime();
        quickSort(array, 0, lines - 1);
        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf", (end - start));
        writeToFile(array, lines, argv[2]);
        
    }

    if(!strcmp(argv[2], "PQ")) {
        printf("Execução com Parallel QuickSort\n");
        start = omp_get_wtime();

        #pragma omp parallel 
        {
            #pragma omp single
            parallel_quickSort(array, 0, lines -1);
        }

        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf", (end - start));
        writeToFile(array, lines, argv[2]);

    }

    if(!strcmp(argv[2], "S")) {
        printf("Execução com Selection Sort\n");
        start = omp_get_wtime();
        selectionSort(array, lines);
        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf", (end - start));
        writeToFile(array, lines, argv[2]);
    }

    if(!strcmp(argv[2], "PS")) {
        printf("Execução com Parallel Selection Sort\n");
        start = omp_get_wtime();
        parallel_selectionSort(array, lines);
        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf", (end - start));
        writeToFile(array, lines, argv[2]);
    }

    return 0;
}


