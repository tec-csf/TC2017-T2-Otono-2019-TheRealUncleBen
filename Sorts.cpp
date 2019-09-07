#include <iostream>
#include <vector> 
#include <algorithm> 
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <bits/stdc++.h>
#define TAM 1000

using namespace std;

void swap(int *xp, int *yp);
void randomize(int a[]);

void bubbleSort(int arr[]);
void CocktailSort(int a[]);
void insertionSort(int arr[]);
void bucketSort(int arr[]);
void countSort(int arr[], int exp);
void selectionSort(int arr[]);
int shellSort(int arr[]);

void heapSort(int arr[]);
void heapify(int arr[], int n, int i);

void merge(int arr[], int l, int m, int r);
void mergeSort(int arr[], int l, int r);

//void treeSort(int arr[]);
//void storeSorted(Node *root, int arr[], int &i);
//Node* insert(Node* node, int key);

int getMax(int arr[]);
void radixsort(int arr[]);

void quickSort(int arr[], int low, int high);
int partition (int arr[], int low, int high);

void printArray(int a[]); 


int main(){
    int *p1, algo;
    p1 = (int *) malloc(TAM*sizeof(int));
	algo = 1;

	randomize(p1);

	auto start = std::chrono::high_resolution_clock::now();
	bubbleSort(p1);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	cout << "La duracion del algoritmo " << algo++ << " fue: " << duration.count() << " millisegundos" << endl;

	randomize(p1);

	start = std::chrono::high_resolution_clock::now();
	insertionSort(p1);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	cout << "La duracion del algoritmo " << algo++ << " fue: " << duration.count() << " millisegundos" << endl;

	randomize(p1);

	start = std::chrono::high_resolution_clock::now();
    CocktailSort(p1);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	cout << "La duracion del algoritmo " << algo++ << " fue: " << duration.count() << " millisegundos" << endl;
	
	/*

	mergeSort(p1, 0, TAM - 1);
	treeSort(p1); 
	radixsort(p1); 
	shellSort(p1);
	selectionSort(p1);
	heapSort(p1);
	quickSort(p1, 0, TAM - 1); 
	
	*/

    free(p1);
    return 0;
}

void printArray(int a[]) { 
	for (int i = 0; i < TAM; i++) 
		printf("%d ", a[i]); 
	printf("\n"); 
} 

void randomize(int a[]){
	srand(1);
    for(int i = 0; i < TAM; i++)
    	a[i] = rand() % (TAM*10);
}
void swap(int *xp, int *yp) { 
	int temp = *xp; 
	*xp = *yp; 
	*yp = temp; 
} 

// A function to implement bubble sort 
void bubbleSort(int arr[]) { 
	int i, j; 
	for (i = 0; i < TAM-1; i++)	 
	
	// Last i elements are already in place 
	for (j = 0; j < TAM-i-1; j++) 
		if (arr[j] > arr[j+1]) 
			swap(&arr[j], &arr[j+1]); 
} 

void CocktailSort(int a[]) { 
	bool swapped = true; 
	int start = 0; 
	int end = TAM - 1; 

	while (swapped) { 
		// reset the swapped flag on entering 
		// the loop, because it might be true from 
		// a previous iteration. 
		swapped = false; 

		// loop from left to right same as 
		// the bubble sort 
		for (int i = start; i < end; ++i) { 
			if (a[i] > a[i + 1]) { 
				swap(a[i], a[i + 1]); 
				swapped = true; 
			} 
		} 

		// if nothing moved, then array is sorted. 
		if (!swapped) 
			break; 

		// otherwise, reset the swapped flag so that it 
		// can be used in the next stage 
		swapped = false; 

		// move the end point back by one, because 
		// item at the end is in its rightful spot 
		--end; 

		// from right to left, doing the 
		// same comparison as in the previous stage 
		for (int i = end - 1; i >= start; --i) { 
			if (a[i] > a[i + 1]) { 
				swap(a[i], a[i + 1]); 
				swapped = true; 
			} 
		} 

		// increase the starting point, because 
		// the last stage would have moved the next 
		// smallest number to its rightful spot. 
		++start; 
	} 
} 

void insertionSort(int arr[]) { 
	int i, key, j; 
	for (i = 1; i < TAM; i++) 
	{ 
		key = arr[i]; 
		j = i - 1; 

		/* Move elements of arr[0..i-1], that are 
		greater than key, to one position ahead 
		of their current position */
		while (j >= 0 && arr[j] > key) { 
			arr[j + 1] = arr[j]; 
			j = j - 1; 
		} 
		arr[j + 1] = key; 
	} 
} 

void bucketSort(int arr[]) { 
	// 1) Create n empty buckets 
	vector<int> b[TAM]; 
	// 2) Put array elements in different buckets 
	cout << "esto: ";
	for (int i=0; i<TAM; i++) { 
		int bi = TAM*arr[i]; // Index in bucket 
		b[bi].push_back(arr[i]); 
	} 
	cout << "jala" << endl;
	// 3) Sort individual buckets 
	for (int i=0; i<TAM; i++) 
	sort(b[i].begin(), b[i].end()); 

	// 4) Concatenate all buckets into arr[] 
	int index = 0; 
	for (int i = 0; i < TAM; i++) {
		for (int j = 0; j < b[i].size(); j++) 
			arr[index++] = b[i][j];
	} 
} 

void Counting_Sort(int A[],int B[]) {
	int k=0;
	for(int i=1;i<=TAM;i++)        
	{
		cin>>A[i];
		if(A[i]>k)
		{
			/*It will modify k if an element 
			occurs whose value is greater than k*/
			k=A[i];              
		}
	}
	int C[k];
	for(int i=0;i<k+1;i++)
	{
		/*It will initialize the C with zero*/
		C[i]=0;
	}
	for(int j=1;j<=TAM;j++)
	{
		/*It will count the occurence of every element x in A 
		and increment it at position x in C*/
		C[A[j]]++;			    
	}
	for(int i=1;i<=k;i++)
	{
		/*It will store the last 
		occurence of the element i */
		C[i]+=C[i-1];            
	}
	for(int j=TAM;j>=1;j--)
	{
		/*It will place the elements at their 
		respective index*/
		B[C[A[j]]]=A[j];          
		/*It will help if an element occurs 
		more than one time*/
		C[A[j]]=C[A[j]]-1;		  
	}
}

void merge(int arr[], int l, int m, int r) { 
	int i, j, k; 
	int n1 = m - l + 1; 
	int n2 = r - m; 

	/* create temp arrays */
	int L[n1], R[n2]; 

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++) 
		L[i] = arr[l + i]; 
	for (j = 0; j < n2; j++) 
		R[j] = arr[m + 1+ j]; 

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray 
	j = 0; // Initial index of second subarray 
	k = l; // Initial index of merged subarray 
	while (i < n1 && j < n2) 
	{ 
		if (L[i] <= R[j]) 
		{ 
			arr[k] = L[i]; 
			i++; 
		} 
		else
		{ 
			arr[k] = R[j]; 
			j++; 
		} 
		k++; 
	} 

	/* Copy the remaining elements of L[], if there 
	are any */
	while (i < n1) 
	{ 
		arr[k] = L[i]; 
		i++; 
		k++; 
	} 

	/* Copy the remaining elements of R[], if there 
	are any */
	while (j < n2) 
	{ 
		arr[k] = R[j]; 
		j++; 
		k++; 
	} 
} 
void mergeSort(int arr[], int l, int r) { 
	if (l < r) 
	{ 
		// Same as (l+r)/2, but avoids overflow for 
		// large l and h 
		int m = l+(r-l)/2; 

		// Sort first and second halves 
		mergeSort(arr, l, m); 
		mergeSort(arr, m+1, r); 

		merge(arr, l, m, r); 
	} 
} 
/*
struct Node 
{ 
	int key; 
	struct Node *left, *right; 
}; 

struct Node *newNode(int item) 
{ 
	struct Node *temp = new Node; 
	temp->key = item; 
	temp->left = temp->right = NULL; 
	return temp; 
} 
void storeSorted(Node *root, int arr[], int &i) 
{ 
	if (root != NULL) 
	{ 
		storeSorted(root->left, arr, i); 
		arr[i++] = root->key; 
		storeSorted(root->right, arr, i); 
	} 
} 
Node* insert(Node* node, int key) 
{ 
	if (node == NULL) return newNode(key); 

	
	if (key < node->key) 
		node->left = insert(node->left, key); 
	else if (key > node->key) 
		node->right = insert(node->right, key); 

	
	return node; 
} 
// This function sorts arr[0..n-1] using Tree Sort 
void treeSort(int arr[]) { 
	struct Node *root = NULL; 

	// Construct the BST 

	root = insert(root, arr[0]); 
	for (int i=1; i<TAM; i++) 
		insert(root, arr[i]); 

	int i = 0; 
	storeSorted(root, arr, i); 
} 
*/
// A utility function to get maximum value in arr[] 
int getMax(int arr[]) 
{ 
	int mx = arr[0]; 
	for (int i = 1; i < TAM; i++) 
		if (arr[i] > mx) 
			mx = arr[i]; 
	return mx; 
} 

void countSort(int arr[], int exp) { 
	int output[TAM]; // output array 
	int i, count[10] = {0}; 

	// Store count of occurrences in count[] 
	for (i = 0; i < TAM; i++) 
		count[ (arr[i]/exp)%10 ]++; 

	// Change count[i] so that count[i] now contains actual 
	// position of this digit in output[] 
	for (i = 1; i < 10; i++) 
		count[i] += count[i - 1]; 

	// Build the output array 
	for (i = TAM - 1; i >= 0; i--) 
	{ 
		output[count[ (arr[i]/exp)%10 ] - 1] = arr[i]; 
		count[ (arr[i]/exp)%10 ]--; 
	} 

	// Copy the output array to arr[], so that arr[] now 
	// contains sorted numbers according to current digit 
	for (i = 0; i < TAM; i++) 
		arr[i] = output[i]; 
} 

// The main function to that sorts arr[] of size n using 
// Radix Sort 
void radixsort(int arr[]) { 
	// Find the maximum number to know number of digits 
	int m = getMax(arr); 

	// Do counting sort for every digit. Note that instead 
	// of passing digit number, exp is passed. exp is 10^i 
	// where i is current digit number 
	for (int exp = 1; m/exp > 0; exp *= 10) 
		countSort(arr,exp); 
} 


int shellSort(int arr[]) { 
	// Start with a big gap, then reduce the gap 
	for (int gap = TAM/2; gap > 0; gap /= 2) 
	{ 
		// Do a gapped insertion sort for this gap size. 
		// The first gap elements a[0..gap-1] are already in gapped order 
		// keep adding one more element until the entire array is 
		// gap sorted 
		for (int i = gap; i < TAM; i += 1) 
		{ 
			// add a[i] to the elements that have been gap sorted 
			// save a[i] in temp and make a hole at position i 
			int temp = arr[i]; 

			// shift earlier gap-sorted elements up until the correct 
			// location for a[i] is found 
			int j;			 
			for (j = i; j >= gap && arr[j - gap] > temp; j -= gap) 
				arr[j] = arr[j - gap]; 
			
			// put temp (the original a[i]) in its correct location 
			arr[j] = temp; 
		} 
	} 
	return 0; 
}

void selectionSort(int arr[]) { 
	int i, j, min_idx; 

	// One by one move boundary of unsorted subarray 
	for (i = 0; i < TAM-1; i++) 
	{ 
		// Find the minimum element in unsorted array 
		min_idx = i; 
		for (j = i+1; j < TAM; j++) 
		if (arr[j] < arr[min_idx]) 
			min_idx = j; 

		// Swap the found minimum element with the first element 
		swap(&arr[min_idx], &arr[i]); 
	} 
} 
 
void heapify(int arr[], int n, int i) { 
	int largest = i; // Initialize largest as root 
	int l = 2*i + 1; // left = 2*i + 1 
	int r = 2*i + 2; // right = 2*i + 2 

	// If left child is larger than root 
	if (l < n && arr[l] > arr[largest]) 
		largest = l; 

	// If right child is larger than largest so far 
	if (r < n && arr[r] > arr[largest]) 
		largest = r; 

	// If largest is not root 
	if (largest != i) 
	{ 
		swap(arr[i], arr[largest]); 

		// Recursively heapify the affected sub-tree 
		heapify(arr, n, largest); 
	} 
} 
void heapSort(int arr[]) { 
	// Build heap (rearrange array) 
	for (int i = TAM / 2 - 1; i >= 0; i--) 
		heapify(arr, TAM, i); 

	// One by one extract an element from heap 
	for (int i=TAM-1; i>=0; i--) 
	{ 
		// Move current root to end 
		swap(arr[0], arr[i]); 

		// call max heapify on the reduced heap 
		heapify(arr, i, 0); 
	} 
} 

int partition (int arr[], int low, int high) { 
	int pivot = arr[high]; // pivot 
	int i = (low - 1); // Index of smaller element 

	for (int j = low; j <= high - 1; j++) 
	{ 
		// If current element is smaller than the pivot 
		if (arr[j] < pivot) 
		{ 
			i++; // increment index of smaller element 
			swap(&arr[i], &arr[j]); 
		} 
	} 
	swap(&arr[i + 1], &arr[high]); 
	return (i + 1); 
} 
void quickSort(int arr[], int low, int high) { 
	if (low < high) 
	{ 
		/* pi is partitioning index, arr[p] is now 
		at right place */
		int pi = partition(arr, low, high); 

		// Separately sort elements before 
		// partition and after partition 
		quickSort(arr, low, pi - 1); 
		quickSort(arr, pi + 1, high); 
	} 
} 
