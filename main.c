#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h> 

#define min 100
struct dataComp { float n; int i; };

#pragma omp declare reduction(compReduction : struct dataComp : omp_out = omp_in.n > omp_out.n ? omp_in : omp_out)

// Escreve array resultante da ordenção no arquivo de output
// O arquivo de output é decidido de acordo com o type
void writeToFile(float* array, int lines, char* type){ 
    if(!strcmp(type, "Q")) {
        FILE* file = fopen("output_q.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    if(!strcmp(type, "PQ")) {
        FILE* file = fopen("output_pq.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    if(!strcmp(type, "S")) {
        FILE* file = fopen("output_s.txt", "w");
        int i;
        for(i = 0; i < lines; i++) fprintf(file, "%f\n", array[i]);
    }

    if(!strcmp(type, "PS")) {
        FILE* file = fopen("output_ps.txt", "w");
        int i;
        for(i = lines - 1; i >= 0; i--) fprintf(file, "%f\n", array[i]);
    }
}

// Conta as linhas de um arquivo
int countLines(char* filePath){ 
    FILE* file = fopen(filePath, "r");
    char c;
    int count = 0;
    for (c = getc(file); c != EOF; c = getc(file))
            if (c == '\n') count++;

    return count;
}

// Lê o conteúdo do arquivo e retorno um array de float com os dados
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

// função genérica para inversão de dois floats em um array
void swap(float *a, float *b) {
    float t = *a;
    *a = *b;
    *b = t;
}

int partition(float array[], int low, int high) {
    // pivot = maior elemento do array
    float pivot = array[high];
    // i do menor elemento 
    // i = indica elemento a esquerda do pivot 
    int i = (low - 1);
    //compara todos elementos do array com o pivot
    for (int j = low; j < high; j++) {
        // Se o elemento for menor que o pivot
        if (array[j] <= pivot) {
            i++;
            // inverte elementos 
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


void selectionSort(float array[], int n)
{
    int i, j;
    int minimo;
 
    //iteracao na parte nao ordenada do array
    for (i = 0; i < n-1; i++)
    {
        minimo = i;
        // itera sobre o array, invertendo os elementos quando necessário 
        for (j = i+1; j < n; j++) if (array[j] < array[minimo]) minimo = j;

        //Inverte menor nor com o primeiro elemento
        if(minimo != i) swap(&array[minimo], &array[i]);
    }
}


void parallel_selectionSort(float* array, int lines){	
    int inicial;
    //itera sobre todas posições do array 
	for(inicial = 0; inicial < lines; inicial++){
        //a cada iteração, define o máximo e seu indice
		struct dataComp data;
        data.n = array[inicial];
        data.i = inicial;

        #pragma omp parallel for reduction(compReduction:data)
		for(int i=inicial +1; i< lines; ++i){
			if(array[i] > data.n){
				data.n = array[i];
				data.i = i;
			}
		}

        // inverte elemento atual com o maior elemento encontrado
		swap(&array[inicial], &array[data.i]);
	}
}


void parallel_quickSort(float * array, int low, int high)
{
    int p;

    if(low < high){ 
        p = partition(array, low, high); 
        #pragma omp task shared(array) if(high - low > min) 
        parallel_quickSort(array, low, p - 1); 
        #pragma omp task shared(array) if(high - low > min)
        parallel_quickSort(array, p + 1, high); 
    }
}


int main (int argc, char** argv){ 

    double start, end;
    int lines = countLines(argv[1]);
    float* array = readFromFile(argv[1], lines);

    omp_set_num_threads(4);

    if(!strcmp(argv[2], "Q")) {
        printf("Execução com QuickSort\n");
        start = omp_get_wtime();
        quickSort(array, 0, lines - 1);
        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf\n", (end - start));
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
        printf("Tempo de Execução: %.5lf\n", (end - start));
        writeToFile(array, lines, argv[2]);

    }

    if(!strcmp(argv[2], "S")) {
        printf("Execução com Selection Sort\n");
        start = omp_get_wtime();
        selectionSort(array, lines);
        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf\n", (end - start));
        writeToFile(array, lines, argv[2]);
    }

    if(!strcmp(argv[2], "PS")) {
        printf("Execução com Parallel Selection Sort\n");
        start = omp_get_wtime();
        parallel_selectionSort(array, lines);
        end = omp_get_wtime();
        printf("Tempo de Execução: %.5lf\n", (end - start));
        writeToFile(array, lines, argv[2]);
    }

    return 0;
}


