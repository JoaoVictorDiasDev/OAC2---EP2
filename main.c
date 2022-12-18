#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printArray(float* array, int lines, int type){ 
    // quick
    if(type == 0) {
        FILE* file = fopen("outputquick.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    FILE* file = fopen("outputselect.txt", "w");
    int i;
    for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);

}

void swap(float *a, float *b) {
    float t = *a;
    *a = *b;
    *b = t;
}

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

int countLines(char* filePath){ 
    FILE* file = fopen(filePath, "r");
    char c;
    int count = 0;
    for (c = getc(file); c != EOF; c = getc(file))
            if (c == '\n') count++;

    return count;
}




int main (int argc, char** argv){ 

    int lines = countLines(argv[1]);
    clock_t start, end;

    //QUICK SORT
    float* arrayQuick = readFromFile(argv[1], lines);
    //------------------------------------------------------------
    double time_quick;
    start = clock();
    quickSort(arrayQuick, 0, lines - 1);
    end = clock();
    time_quick = ((double) (end - start)) / CLOCKS_PER_SEC; 
    //-------------------------------------------------------------
    printArray (arrayQuick, lines, 0);


    //SELECTION SORT
    float* arraySelect = readFromFile(argv[1], lines);
    //------------------------------------------------------------
    double time_select;
    start = clock();
    selectionSort(arraySelect, lines);
    end = clock();
    time_select = ((double) (end - start)) / CLOCKS_PER_SEC; 
    //------------------------------------------------------------
    printArray (arraySelect, lines, 1);

    printf("Temp de execucao:\nQuickSort -> %lf\nSelectionSort -> %lf\n", time_quick, time_select);

    return 0;
}


