#include<iostream>
#include<pthread.h>
#include<cmath>
#include<sys/time.h>
#include<string.h>
using namespace std;

// Structure for passing arguments for parallel threads
struct quicksort_structure {
  int *a;
  int a_start;
  int low;
  int high;
  int pivot;
  int threadId;
  int noOfProcessorsInGroup;
  int maxThreads;
  pthread_barrier_t *barrier;
  int *less_sum;
  int *high_sum;
  int *cumulative_less_sum;
  int *cumulative_high_sum;
  int *temp_array;
};

// Structure for recursive calls of the main parallel sort method
struct main_arg_structure {
    int *a;
    int low;
    int high;
    int size;
    int noOfThreads;
};

// Method to swap two elements of an array
void swap(int *a, int *b) {
    int temp = *b;
    *b = *a;
    *a = temp;
}

// Method for quick sort increasing order sort comparator
int quick_sort_comparator (const void * c, const void * d) {
    int *temp1 = (int *) c;
    int *temp2 = (int *) d;
    return (*temp1 - *temp2);
}

// Method that takes a pivot and finds its position in the sub array and calculates number of elements less than and greater than that
void localArrange(int *a, int low, int high, int pivot, int threadId, int *less_sum, int *high_sum, int *cumulative_less_sum, int *cumulative_high_sum) {
    int k = high;
    for (int i = high; i >= low; i--) {
        if (a[i] >= pivot) {
            swap(&a[i], &a[k]);
            k--;
        }
    }

    less_sum[threadId] = k - low + 1;
    high_sum[threadId] = high - k;
    cumulative_less_sum[threadId] = less_sum[threadId];
    cumulative_high_sum[threadId] = high_sum[threadId];
}

// Method that takes cumulative prefix sum, and arranges the elements of the subarray in the global array based on the common pivot
void globalArrange(int *a, int low, int high, int threadId, int noOfThreads, int *cumulative_less_sum, int *cumulative_high_sum, int *less_sum,
int *temp_array, pthread_barrier_t *local_barrier, int a_start) {
    int local_low_left = cumulative_less_sum[threadId];
    int local_high_left = cumulative_less_sum[threadId+1];

    int local_low_right = cumulative_less_sum[noOfThreads] + cumulative_high_sum[threadId];
    int local_high_right = cumulative_less_sum[noOfThreads] + cumulative_high_sum[threadId+1];

    int total = local_high_left - local_low_left;

    int j = 0;
    for(int i = local_low_left; i < local_high_left; i++) {
        temp_array[i] = a[low + j++];
    }

    j = 0;
    for(int i = local_low_right; i < local_high_right; i++) {
        temp_array[i] = a[low + less_sum[threadId] + j++];
    }

    pthread_barrier_wait(local_barrier);

    int size = high-low+1;

    for(int i=0;i<size;i++) {
        a[low+i] = temp_array[low+i-a_start];
    }

}

// Up pass to calculate cumulative prefix sum
void tree_up_pass(int *cumulative_sum, int level, int myStart) {
    int firstIndex = myStart + int(pow(2,level)) - 1;
    int secondIndex = myStart + int(pow(2,(level+1))) - 1;
    cumulative_sum[secondIndex] = cumulative_sum[firstIndex] + cumulative_sum[secondIndex];
    return;
}

// Down pass to calculate cumulative prefix sum
void tree_down_pass(int *cumulative_sum, int level, int myStart) {
    int firstIndex = myStart + int(pow(2,level)) - 1;
    int secondIndex = myStart + int(pow(2,level+1)) - 1;

    int temp_value = cumulative_sum[firstIndex];
    cumulative_sum[firstIndex] = cumulative_sum[secondIndex];
    cumulative_sum[secondIndex] = temp_value + cumulative_sum[secondIndex];
    return;
}

// Method that performs the local partition and global rearrangement for individual threads
void *parallel_quicksort(void *quicksort_arg) {
    quicksort_structure *local_arg = (quicksort_structure *) quicksort_arg; 
    int *a = local_arg->a;
    int a_start = local_arg->a_start;
    int low = local_arg->low;
    int high = local_arg->high;
    int pivot = local_arg->pivot;
    int threadId = local_arg->threadId;
    int noOfThreads = local_arg->noOfProcessorsInGroup;
    int maxThreads = local_arg->maxThreads;
    pthread_barrier_t *local_barrier = local_arg->barrier;
    int *less_sum = local_arg->less_sum;
    int *high_sum = local_arg->high_sum;
    int *temp_array = local_arg->temp_array;
    int *cumulative_less_sum = local_arg->cumulative_less_sum;
    int *cumulative_high_sum = local_arg->cumulative_high_sum;

    localArrange(a, low, high, pivot, threadId, less_sum, high_sum, cumulative_less_sum, cumulative_high_sum);

    pthread_barrier_wait(local_barrier);

    // Compute cumulative prefix sum using the tree structure in log (P) where P is the number of threads
    for (int level = 0; level <= log2(maxThreads)-1; level++) {
        int myStart = threadId * int(pow(2,(level+1)));
        if(myStart < maxThreads) {
            tree_up_pass(cumulative_less_sum, level, myStart);
            tree_up_pass(cumulative_high_sum, level, myStart);
        }
        pthread_barrier_wait(local_barrier);
    }

    if(threadId == 0) {
        cumulative_less_sum[maxThreads-1] = 0;
        cumulative_high_sum[maxThreads-1] = 0;
    }
    
    pthread_barrier_wait(local_barrier);
    
    for (int level = log2(maxThreads)-1; level >= 0; level--) {
        int myStart = threadId * int(pow(2,level+1));
        if(myStart < maxThreads) {
            tree_down_pass(cumulative_less_sum, level, myStart);
            tree_down_pass(cumulative_high_sum, level, myStart);
        }
        pthread_barrier_wait(local_barrier);
    }

    if(threadId == 0) {
        if(noOfThreads == maxThreads) {
            cumulative_less_sum[maxThreads] = cumulative_less_sum[maxThreads-1] + less_sum[maxThreads-1];
            cumulative_high_sum[maxThreads] = cumulative_high_sum[maxThreads-1] + high_sum[maxThreads-1];
        } else {
            cumulative_less_sum[maxThreads] = cumulative_less_sum[maxThreads-1];
            cumulative_high_sum[maxThreads] = cumulative_high_sum[maxThreads-1];
        }
    }

    pthread_barrier_wait(local_barrier);

    globalArrange(a, low, high, threadId, noOfThreads, cumulative_less_sum, cumulative_high_sum, less_sum, temp_array, local_barrier, a_start);

    pthread_exit(NULL);
    return NULL;
}

// Method to create the main structure argument
main_arg_structure* create_main_arg_struct(int *a, int low, int high, int noOfThreads, int size) {
    main_arg_structure *new_arg = new main_arg_structure;
    new_arg->a = a;
    new_arg->low = low;
    new_arg->high = high;
    new_arg->noOfThreads = noOfThreads;
    new_arg->size = size;
    return new_arg;
}

// Method to create the quick sort structure argument
quicksort_structure* create_quick_sort_struct(int *a, int a_start, int start_index, int end_index, int pivot, int threadId, int noOfThreads,
int maxThreads, pthread_barrier_t *barrier, int *less_sum, int *high_sum, int *temp_array, int *cumulative_less_sum,
int *cumulative_high_sum) {
    quicksort_structure* new_arg = new quicksort_structure;
    new_arg->a = a;
    new_arg->a_start = a_start;
    new_arg->low = start_index;
    new_arg->high = end_index;
    new_arg->pivot = pivot;
    new_arg->threadId = threadId;
    new_arg->noOfProcessorsInGroup = noOfThreads;
    new_arg->maxThreads = maxThreads;
    new_arg->barrier = barrier;
    new_arg->less_sum = less_sum;
    new_arg->high_sum = high_sum;
    new_arg->temp_array = temp_array;
    new_arg->cumulative_less_sum = cumulative_less_sum;
    new_arg->cumulative_high_sum = cumulative_high_sum;
    return new_arg;
}

// Method for parallel quick sort
void *quick_sort(void *main_arg) {
    main_arg_structure *local_main_arg = (main_arg_structure *) main_arg;
    int noOfThreads = local_main_arg->noOfThreads;
    int low = local_main_arg->low;
    int high = local_main_arg->high;
    int size = local_main_arg->size;
    int *a = local_main_arg->a;

    if(low < high) {

        if(noOfThreads <= 1) {
            qsort(&a[low], high - low + 1, sizeof(int), quick_sort_comparator);
            pthread_exit(NULL);
            return NULL;
        }

        pthread_t worker_threads[noOfThreads];
        
        pthread_barrier_t barrier;
        pthread_barrier_init(&barrier, NULL, noOfThreads);

        int pivot = a[low + rand() % size];

        int maxThreads = int(pow(2, ceil(log2(noOfThreads))));

        int *less_sum = new int[maxThreads];
        int *high_sum = new int[maxThreads];
        int *temp_array = new int[size];
        int *cumulative_less_sum = new int[maxThreads+1];
        int *cumulative_high_sum = new int[maxThreads+1];
        int block = size / noOfThreads;
        int remainder = size % noOfThreads;
        int start_index = low;

        for(int i = 0; i<noOfThreads; i++) {
            int size = block + (i < remainder ? 1 : 0);
            int end_index = start_index + size - 1;
            quicksort_structure* parallel_quicksort_arg = create_quick_sort_struct(a, low, start_index, end_index, pivot, i, noOfThreads,
            maxThreads, &barrier, less_sum, high_sum, temp_array, cumulative_less_sum, cumulative_high_sum);
            pthread_create(&worker_threads[i], NULL, parallel_quicksort, parallel_quicksort_arg);
            start_index = start_index + size;
        }

        for(int i=0; i<noOfThreads; i++) {
            pthread_join(worker_threads[i], NULL);
        }

        pthread_barrier_destroy(&barrier);

        int pivot_index = cumulative_less_sum[noOfThreads] + low;

        // Compute the number of threads of each sub array
        int left_threads = ceil(((cumulative_less_sum[noOfThreads] * double(noOfThreads)) / size ) + 0.5);
        left_threads = (left_threads > noOfThreads) ? noOfThreads : left_threads;
        int right_threads = noOfThreads - left_threads;

        int left_size = pivot_index - low;
        int right_size = high - pivot_index + 1;

        if( (left_size < right_size && left_threads > right_threads) 
         || (right_size < left_size && right_threads > left_threads)) {
            swap(&left_threads, &right_threads);
        }

        main_arg_structure* left_arg = create_main_arg_struct(a, low, pivot_index-1, left_threads, left_size);
        main_arg_structure* right_arg = create_main_arg_struct(a, pivot_index, high, right_threads, right_size);

        delete[] less_sum;
        delete[] high_sum;
        delete[] cumulative_less_sum;
        delete[] cumulative_high_sum;
        delete[] temp_array;
        
        // Sort left and right subarrays parallelly
        pthread_t left_thread;
        pthread_t right_thread;
        pthread_create(&left_thread, NULL, quick_sort, left_arg);
        pthread_create(&right_thread, NULL, quick_sort, right_arg);
        pthread_join(left_thread, NULL);
        pthread_join(right_thread, NULL);
    }

    pthread_exit(NULL);
    return NULL;

}

void performSort(int noOfThreads) {
    int array_sizes[] = {100, 1000, 10000, 100000, 1000000, 10000000};

    for(int i = 0; i < sizeof(array_sizes)/sizeof(int); i++) {
        struct timeval start, end;  
        int n = array_sizes[i];
        int *a = new int[n];
        int *temp = new int[n];

        for(int i = 0; i<n;i++) {
            a[i] = rand() % n;
        }

        memcpy(temp, a, n*sizeof(int));

        gettimeofday(&start, NULL);
        qsort(temp, n, sizeof(int), quick_sort_comparator);
        gettimeofday(&end, NULL);

        long startTime = (start.tv_sec * 1000000) + start.tv_usec;
        long endTime = (end.tv_sec * 1000000) + end.tv_usec;

        double serialTime = double(endTime - startTime) / 1000000;

        printf("Time for serial quicksort for size %d: %lf seconds \n", array_sizes[i], serialTime);
        
        main_arg_structure* main_arg = create_main_arg_struct(a, 0, n-1, noOfThreads, n);
        pthread_t main_thread;
        gettimeofday(&start, NULL);
        pthread_create(&main_thread, NULL, quick_sort, main_arg);
        pthread_join(main_thread, NULL);
        gettimeofday(&end, NULL);

        startTime = (start.tv_sec * 1000000) + start.tv_usec;
        endTime = (end.tv_sec * 1000000) + end.tv_usec;

        double parallelTime = double(endTime - startTime) / 1000000;

        printf("Time for parallel quicksort with %d threads for size %d: %lf seconds \n", noOfThreads, array_sizes[i], parallelTime);

        double speedUp = double(serialTime) / double(parallelTime);

        printf("Speedup for %d threads with size %d: %lf \n", noOfThreads, array_sizes[i], speedUp);

        for(int i=0;i<n-1;i++) {
            if(a[i] > a[i+1]) {
                printf("Not sorted\n");
            }
        }

        delete[] a;
        delete[] temp;
    }
}

int main() {

    int num_threads[] = {4, 8, 16, 32, 64, 128};
    srand(time(NULL));

    for(int i = 0; i < sizeof(num_threads)/sizeof(int); i++) {
        struct timeval start, end;  
        int n = 10000000;
        int *a = new int[n];
        int *temp = new int[n];

        printf("Timings for array size %d \n", n);

        for(int i = 0; i<n;i++) {
            a[i] = rand() % n;
        }

        memcpy(temp, a, n*sizeof(int));

        gettimeofday(&start, NULL);
        qsort(temp, n, sizeof(int), quick_sort_comparator);
        gettimeofday(&end, NULL);

        long startTime = (start.tv_sec * 1000000) + start.tv_usec;
        long endTime = (end.tv_sec * 1000000) + end.tv_usec;

        double serialTime = double(endTime - startTime) / 1000000;

        printf("Time for serial quicksort: %lf seconds \n", serialTime);

        main_arg_structure* main_arg = create_main_arg_struct(a, 0, n-1, num_threads[i], n);
        pthread_t main_thread;
        gettimeofday(&start, NULL);
        pthread_create(&main_thread, NULL, quick_sort, main_arg);
        pthread_join(main_thread, NULL);
        gettimeofday(&end, NULL);

        startTime = (start.tv_sec * 1000000) + start.tv_usec;
        endTime = (end.tv_sec * 1000000) + end.tv_usec;

        double parallelTime = double(endTime - startTime) / 1000000;

        printf("Time for parallel quicksort with %d threads: %lf seconds \n", num_threads[i], parallelTime);

        double speedUp = double(serialTime) / double(parallelTime);

        printf("Speedup for %d threads: %lf \n", num_threads[i], speedUp);

        for(int i=0;i<n-1;i++) {
            if(a[i] > a[i+1]) {
                printf("Not sorted\n");
            }
        }

        delete[] a;
        delete[] temp;
    }

    printf("\n");

    printf("############################################ \n");

    printf("\n");

    performSort(16);

    printf("\n");

    printf("############################################ \n");

    printf("\n");

    performSort(32);

    return 1;
}